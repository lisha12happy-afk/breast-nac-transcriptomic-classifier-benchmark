param(
    [string]$OutDir = "manifests"
)

$ErrorActionPreference = "Stop"

function Normalize-Token {
    param([string]$Value)
    if ($null -eq $Value) { return "" }
    $v = [regex]::Replace(([string]$Value).Trim(), '\s+', ' ')
    $v = $v.Trim('"')
    return $v.Trim()
}

function Get-SeriesPrefix {
    param([string]$Accession)
    if ($Accession.Length -le 3) { return "${Accession}nnn" }
    return $Accession.Substring(0, $Accession.Length - 3) + "nnn"
}

function Get-GEO-MatrixFiles {
    param([string]$Accession)
    $prefix = Get-SeriesPrefix $Accession
    $url = "https://ftp.ncbi.nlm.nih.gov/geo/series/$prefix/$Accession/matrix/"
    $r = Invoke-WebRequest -Uri $url -UseBasicParsing -TimeoutSec 60
    $files = @()
    foreach ($link in $r.Links) {
        $href = [string]$link.href
        if ($href -match '\.txt\.gz$') {
            $files += ($href -split '/')[-1]
        }
    }
    return @($files | Sort-Object -Unique)
}

function Read-GEO-MatrixHeader {
    param(
        [string]$Accession,
        [string]$MatrixFile
    )
    $prefix = Get-SeriesPrefix $Accession
    $url = "https://ftp.ncbi.nlm.nih.gov/geo/series/$prefix/$Accession/matrix/$MatrixFile"
    $req = [System.Net.WebRequest]::Create($url)
    $req.Timeout = 90000
    $resp = $req.GetResponse()
    try {
        $stream = $resp.GetResponseStream()
        $gz = New-Object IO.Compression.GzipStream($stream, [IO.Compression.CompressionMode]::Decompress)
        $sr = New-Object IO.StreamReader($gz)
        $lines = New-Object System.Collections.Generic.List[string]
        while (($line = $sr.ReadLine()) -ne $null) {
            if ($line -like "!Sample_*") { $lines.Add($line) }
            if ($line -eq "!series_matrix_table_begin") { break }
        }
        $sr.Close()
        $gz.Close()
        $stream.Close()
        return @{ url = $url; lines = @($lines) }
    }
    finally {
        $resp.Close()
    }
}

function Parse-GEO-Samples {
    param(
        [string]$Accession,
        [string]$MatrixFile,
        [string[]]$Lines
    )
    $sampleCount = 0
    foreach ($line in $Lines) {
        if (-not ($line -like "!Sample_*")) { continue }
        $parts = $line -split "`t"
        if ($parts.Count -gt 1 -and ($parts.Count - 1) -gt $sampleCount) {
            $sampleCount = $parts.Count - 1
        }
    }

    $sampleMetas = New-Object System.Collections.Generic.List[object]
    for ($i = 0; $i -lt $sampleCount; $i++) {
        $sampleMetas.Add([ordered]@{}) | Out-Null
    }

    foreach ($line in $Lines) {
        if (-not ($line -like "!Sample_*")) { continue }
        $parts = $line -split "`t"
        if ($parts.Count -lt 2) { continue }
        $field = $parts[0].TrimStart("!")
        $values = @($parts | Select-Object -Skip 1 | ForEach-Object { Normalize-Token $_ })

        for ($i = 0; $i -lt $sampleCount; $i++) {
            if ($i -ge $values.Count) { continue }
            $value = Normalize-Token $values[$i]
            if ($value -eq "") { continue }
            $meta = [System.Collections.IDictionary]$sampleMetas.Item($i)

            if ($field -like "Sample_characteristics_ch*") {
                $key = ""
                $val = ""
                if ($value -match '^([^:]+)\s*:\s*(.*)$') {
                    $key = Normalize-Token $matches[1]
                    $val = Normalize-Token $matches[2]
                }
                elseif ($value -match '^([^;]+);(.+)$') {
                    # Some GEO records use semicolon-delimited characteristics without a colon.
                    $key = Normalize-Token $matches[1]
                    $val = Normalize-Token $matches[2]
                }
                if ($key -ne "") {
                    if ($meta.Contains($key)) {
                        $old = Normalize-Token $meta[$key]
                        if ($val -ne "" -and $old -ne $val -and $old -notmatch [regex]::Escape($val)) {
                            $meta[$key] = (($old, $val | Where-Object { $_ -ne "" } | Select-Object -Unique) -join " | ")
                        }
                    }
                    else {
                        $meta[$key] = $val
                    }
                }
            }
            else {
                if ($meta.Contains($field)) {
                    $old = Normalize-Token $meta[$field]
                    if ($old -ne $value -and $old -notmatch [regex]::Escape($value)) {
                        $meta[$field] = (($old, $value | Where-Object { $_ -ne "" } | Select-Object -Unique) -join " | ")
                    }
                }
                else {
                    $meta[$field] = $value
                }
            }
        }
    }

    $sampleRows = @()
    for ($i = 0; $i -lt $sampleCount; $i++) {
        $sampleRows += [pscustomobject]@{
            accession = $Accession
            matrix_file = $MatrixFile
            meta = $sampleMetas.Item($i)
        }
    }
    return $sampleRows
}

function First-Meta {
    param(
        [System.Collections.IDictionary]$Meta,
        [string[]]$Keys
    )
    foreach ($k in $Keys) {
        $hasKey = $false
        if ($Meta -is [System.Collections.IDictionary]) { $hasKey = $Meta.Contains($k) }
        if (-not $hasKey -and ($Meta.PSObject.Methods.Name -contains "ContainsKey")) { $hasKey = $Meta.ContainsKey($k) }
        if ($hasKey) {
            $v = Normalize-Token $Meta[$k]
            if ($v -ne "" -and $v -notin @("NA", "na", "N/A", "unknown", "Unknown")) { return $v }
        }
    }
    return ""
}

function Map-Endpoint {
    param(
        [string]$Dataset,
        [System.Collections.IDictionary]$Meta
    )
    $rawField = ""
    $raw = ""
    switch ($Dataset) {
        "GSE25066" { $rawField = "pathologic_response_pcr_rd"; $raw = First-Meta $Meta @($rawField) }
        "GSE20271" { $rawField = "pcr or rd"; $raw = First-Meta $Meta @($rawField) }
        "GSE32646" { $rawField = "pathologic response pcr ncr"; $raw = First-Meta $Meta @($rawField) }
        "GSE41998" { $rawField = "pcr"; $raw = First-Meta $Meta @($rawField) }
        "GSE50948" { $rawField = "pcr"; $raw = First-Meta $Meta @($rawField) }
        "GSE66305" { $rawField = "pcr (1=yes)"; $raw = First-Meta $Meta @($rawField) }
        "GSE163882" { $rawField = "response to nac"; $raw = First-Meta $Meta @($rawField) }
        "GSE194040" { $rawField = "pcr"; $raw = First-Meta $Meta @($rawField) }
        "GSE106977" { $rawField = "pathological complete response"; $raw = First-Meta $Meta @($rawField) }
        "GSE109710" { $rawField = "pcr"; $raw = First-Meta $Meta @($rawField) }
        "GSE130786" { $rawField = "drug response"; $raw = First-Meta $Meta @($rawField) }
        default { $rawField = ""; $raw = "" }
    }
    $norm = ""
    $bin = ""
    $v = $raw.ToLowerInvariant()
    if ($Dataset -eq "GSE41998") {
        if ($v -eq "yes") {
            return [pscustomobject]@{ raw_field = $rawField; raw_value = $raw; endpoint = "pCR"; binary = "1" }
        }
        elseif ($v -eq "no") {
            return [pscustomobject]@{ raw_field = $rawField; raw_value = $raw; endpoint = "RD"; binary = "0" }
        }
        else {
            return [pscustomobject]@{ raw_field = $rawField; raw_value = $raw; endpoint = ""; binary = "" }
        }
    }
    if ($v -in @("pcr", "pcr ", "yes", "1", "pcr=yes", "pcr=1", "cr", "pcr+npcr mamma")) {
        $norm = "pCR"; $bin = "1"
    }
    elseif ($v -in @("rd", "ncr", "no", "0", "non-pcr", "nonpcr", "npcr", "residual disease", "no pcr")) {
        $norm = "RD"; $bin = "0"
    }
    elseif ($v -match '^pcr$|^yes$|^1$') {
        $norm = "pCR"; $bin = "1"
    }
    elseif ($v -match '^rd$|^ncr$|^no$|^0$|npcr|residual') {
        $norm = "RD"; $bin = "0"
    }
    return [pscustomobject]@{ raw_field = $rawField; raw_value = $raw; endpoint = $norm; binary = $bin }
}

function Map-Status {
    param(
        [string]$Raw,
        [string]$Default = ""
    )
    $rawNorm = Normalize-Token $Raw
    if ($rawNorm -eq "" -and $Default -ne "") { return $Default }
    $v = $rawNorm.ToLowerInvariant()
    if ($v -in @("p", "positive", "pos", "1", "yes", "y", "+", "ep", "pp")) { return "Positive" }
    if ($v -in @("n", "negative", "neg", "0", "no", "en", "pn", "-")) { return "Negative" }
    if ($v -match 'positive|pos|\+|^p$') { return "Positive" }
    if ($v -match 'negative|neg|\-|^n$') { return "Negative" }
    if ($rawNorm -eq "" -or $v -in @("na", "n/a", "unknown", "indeterminate")) { return "" }
    return $rawNorm
}

function Derive-Subtype {
    param(
        [string]$Dataset,
        [string]$ER,
        [string]$PR,
        [string]$HER2,
        [string]$HR
    )
    if ($Dataset -eq "GSE106977") { return "TNBC_by_study" }
    $hrVal = $HR
    if ($hrVal -eq "") {
        if ($ER -eq "Positive" -or $PR -eq "Positive") { $hrVal = "Positive" }
        elseif ($ER -eq "Negative" -and $PR -eq "Negative") { $hrVal = "Negative" }
    }
    if ($HER2 -eq "Positive") { return "HER2_positive" }
    if ($HER2 -eq "Negative" -and $hrVal -eq "Negative") { return "TNBC_or_ER_PR_HER2_negative" }
    if ($HER2 -eq "Negative" -and $hrVal -eq "Positive") { return "HR_positive_HER2_negative" }
    if ($hrVal -eq "Negative") { return "HR_negative_HER2_unknown" }
    if ($hrVal -eq "Positive") { return "HR_positive_HER2_unknown" }
    return "unknown"
}

function Get-StatusFields {
    param(
        [string]$Dataset,
        [System.Collections.IDictionary]$Meta
    )
    $studyHER2 = ""
    $studyER = ""
    $studyPR = ""
    $studyHR = ""
    if ($Dataset -in @("GSE50948","GSE66305","GSE109710","GSE130786")) { $studyHER2 = "Positive" }
    if ($Dataset -eq "GSE106977") { $studyER = "Negative"; $studyPR = "Negative"; $studyHER2 = "Negative"; $studyHR = "Negative" }

    $erRaw = First-Meta $Meta @("er_status_ihc", "er status", "er status ihc", "er", "ER", "er.status (1", "er status (ihc staining results)")
    $prRaw = First-Meta $Meta @("pr_status_ihc", "pr status", "pr status ihc", "pr", "PR", "pr.status (1", "pr_status (ihc staining results)", "progesterone receptor status")
    $her2Raw = First-Meta $Meta @("her2_status", "her2 status", "her2 status fish", "her 2 status", "her2", "HER2", "her2stat", "her2 receptor status")
    $hrRaw = First-Meta $Meta @("hr", "hormone receptor status")

    $er = Map-Status $erRaw $studyER
    $pr = Map-Status $prRaw $studyPR
    $her2 = Map-Status $her2Raw $studyHER2
    $hr = Map-Status $hrRaw $studyHR
    if ($hr -eq "" -and ($er -eq "Positive" -or $pr -eq "Positive")) { $hr = "Positive" }
    elseif ($hr -eq "" -and $er -eq "Negative" -and $pr -eq "Negative") { $hr = "Negative" }

    return [pscustomobject]@{
        er_raw = $erRaw; er = $er
        pr_raw = $prRaw; pr = $pr
        her2_raw = $her2Raw; her2 = $her2
        hr_raw = $hrRaw; hr = $hr
        subtype = Derive-Subtype $Dataset $er $pr $her2 $hr
    }
}

function Check-SupplementaryListing {
    param([string]$Accession)
    $prefix = Get-SeriesPrefix $Accession
    $url = "https://ftp.ncbi.nlm.nih.gov/geo/series/$prefix/$Accession/suppl/"
    try {
        $r = Invoke-WebRequest -Uri $url -UseBasicParsing -TimeoutSec 45
        $files = @()
        foreach ($link in $r.Links) {
            $href = [string]$link.href
            if ($href -and $href -notmatch '^\.\.?/?$' -and $href -notmatch 'vulnerability') {
                $files += (($href -split '/')[-1])
            }
        }
        $unique = @($files | Where-Object { $_ } | Sort-Object -Unique)
        $rawHint = if (($unique -join ';') -match '\.CEL|\.fastq|\.fq|\.bam|raw|counts|tpm|fpkm') { "yes_or_possible_supplementary" } elseif ($unique.Count -gt 0) { "supplementary_present_unclear_raw" } else { "no_supplementary_files_listed" }
        return [pscustomobject]@{ url = $url; raw_hint = $rawHint; file_count = $unique.Count; example_files = (($unique | Select-Object -First 5) -join ";") }
    }
    catch {
        return [pscustomobject]@{ url = $url; raw_hint = "not_confirmed"; file_count = ""; example_files = ""; error = $_.Exception.Message }
    }
}

$datasetConfigs = @(
    [pscustomobject]@{ accession="GSE25066"; role="discovery_model_reconstruction_base"; story="full_benchmark"; include_primary_analysis="yes"; overlap_flag="superset_of_GSE25055_GSE25065; do_not_combine_with_GSE25055_or_GSE25065"; notes="chosen over GSE25055 to use the larger Hatzis/MDACC-USO series" },
    [pscustomobject]@{ accession="GSE20271"; role="external_validation_1"; story="full_benchmark"; include_primary_analysis="yes"; overlap_flag="I-SPY1 independent validation; check against I-SPY2 separately"; notes="" },
    [pscustomobject]@{ accession="GSE32646"; role="external_validation_2_and_TNBC_ER_negative"; story="full_benchmark;subtype_TNBC_ER_negative"; include_primary_analysis="yes"; overlap_flag="independent ER-negative cohort"; notes="" },
    [pscustomobject]@{ accession="GSE41998"; role="external_validation_3"; story="full_benchmark"; include_primary_analysis="yes"; overlap_flag="independent trial; treatment arm heterogeneity"; notes="" },
    [pscustomobject]@{ accession="GSE50948"; role="external_validation_4_and_HER2_positive"; story="full_benchmark;subtype_HER2_positive"; include_primary_analysis="yes"; overlap_flag="NOAH HER2-positive trial; independent from CHER-LOB/TRYPHAENA"; notes="" },
    [pscustomobject]@{ accession="GSE66305"; role="external_validation_5_and_HER2_positive"; story="full_benchmark;subtype_HER2_positive"; include_primary_analysis="yes"; overlap_flag="CHER-LOB HER2-positive trial"; notes="" },
    [pscustomobject]@{ accession="GSE163882"; role="external_validation_6_RNAseq_FFPE"; story="full_benchmark"; include_primary_analysis="yes"; overlap_flag="RNA-seq FFPE multi-center; platform stress for microarray-trained models"; notes="" },
    [pscustomobject]@{ accession="GSE194040"; role="large_stress_test_I-SPY2"; story="full_benchmark_stress;subtype_TNBC_ER_negative;subtype_HER2_positive"; include_primary_analysis="yes_stress_test"; overlap_flag="I-SPY2-990 mRNA component; do_not_combine_with_GSE196096_or_I-SPY2_subseries"; notes="two matrix files represent different platforms/batches" },
    [pscustomobject]@{ accession="GSE106977"; role="subtype_validation_TNBC"; story="subtype_TNBC_ER_negative"; include_primary_analysis="subtype_only"; overlap_flag="TNBC-specific subtype validation"; notes="" },
    [pscustomobject]@{ accession="GSE109710"; role="subtype_validation_HER2_positive"; story="subtype_HER2_positive"; include_primary_analysis="subtype_only"; overlap_flag="TRYPHAENA HER2-positive NanoString panel; not whole transcriptome"; notes="use for gene-availability sensitivity and HER2+ subgroup" },
    [pscustomobject]@{ accession="GSE130786"; role="subtype_validation_HER2_positive"; story="subtype_HER2_positive"; include_primary_analysis="subtype_only"; overlap_flag="HER2-positive TCH/TCL/TCHL baseline subset"; notes="use baseline only; do not combine with GSE130787 treatment samples" }
)

New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

$allRows = New-Object System.Collections.Generic.List[object]
$datasetRows = New-Object System.Collections.Generic.List[object]
$dictionaryRows = New-Object System.Collections.Generic.List[object]
$errors = New-Object System.Collections.Generic.List[object]

foreach ($cfg in $datasetConfigs) {
    Write-Host "Processing $($cfg.accession) ..."
    try {
        $matrixFiles = Get-GEO-MatrixFiles $cfg.accession
        $supp = Check-SupplementaryListing $cfg.accession
        foreach ($mf in $matrixFiles) {
            $read = Read-GEO-MatrixHeader $cfg.accession $mf
            $sampleRows = Parse-GEO-Samples $cfg.accession $mf $read.lines
            $platform = ""
            if ($mf -match 'GPL\d+') { $platform = $matches[0] }
            else {
                $gplLine = $read.lines | Where-Object { $_ -like "!Sample_platform_id*" } | Select-Object -First 1
                if ($gplLine) {
                    $vals = @($gplLine -split "`t" | Select-Object -Skip 1 | ForEach-Object { Normalize-Token $_ } | Select-Object -Unique)
                    $platform = ($vals -join ";")
                }
            }

            foreach ($sr in $sampleRows) {
                $m = $sr.meta
                $endpoint = Map-Endpoint $cfg.accession $m
                $status = Get-StatusFields $cfg.accession $m
                $geo = First-Meta $m @("Sample_geo_accession")
                $title = First-Meta $m @("Sample_title")
                $source = First-Meta $m @("source", "Sample_source_name_ch1")
                $patientId = First-Meta $m @("patient id", "patient identifier", "sample id", "Sample_title")
                $timepointRaw = First-Meta $m @("biopsy time", "time point", "timepoint", "treatment time", "visit")
                $timepoint = ""
                if ($cfg.accession -in @("GSE194040","GSE25066","GSE20271","GSE32646","GSE41998","GSE50948","GSE66305","GSE106977","GSE109710","GSE130786","GSE163882")) { $timepoint = "pretreatment_or_baseline_by_series_design" }
                if ($timepointRaw -match 'on|post|week|mid|resection|treat') { $timepoint = $timepointRaw }
                $treatmentArm = First-Meta $m @("arm", "treatment arm", "treatment group", "cher-lob arm", "arm description", "treatment.type", "treatment received (1=fac, 2=t/fac)", "preoperative treatment", "type_taxane")
                $include = $cfg.include_primary_analysis
                $exclusion = ""
                if ($endpoint.binary -eq "") { $include = "no"; $exclusion = "endpoint_missing_or_unmapped" }
                if ($timepoint -match 'post|on-treatment|mid|resection|week|treatment' -and $cfg.accession -notin @("GSE194040")) {
                    # Selected datasets are intended as baseline. Keep this conservative in case a mixed matrix is added later.
                    if ($timepoint -notmatch 'baseline|pretreatment|pre-treatment') { $include = "no"; $exclusion = "not_pretreatment" }
                }
                $allRows.Add([pscustomobject]@{
                    dataset_accession = $cfg.accession
                    matrix_file = $mf
                    gsm_sample_id = $geo
                    sample_title = $title
                    patient_id = $patientId
                    source_or_center = $source
                    timepoint_raw = $timepointRaw
                    timepoint_standard = $timepoint
                    endpoint_raw_field = $endpoint.raw_field
                    endpoint_raw_value = $endpoint.raw_value
                    endpoint_standard = $endpoint.endpoint
                    endpoint_binary_pcr1_rd0 = $endpoint.binary
                    er_raw = $status.er_raw
                    er_status = $status.er
                    pr_raw = $status.pr_raw
                    pr_status = $status.pr
                    her2_raw = $status.her2_raw
                    her2_status = $status.her2
                    hr_raw = $status.hr_raw
                    hr_status = $status.hr
                    subtype_group = $status.subtype
                    treatment_arm = $treatmentArm
                    platform = $platform
                    expression_type = if ($cfg.accession -eq "GSE163882") { "RNA-seq_FFPE" } elseif ($cfg.accession -eq "GSE109710") { "NanoString_targeted_panel" } else { "microarray_or_array" }
                    normalized_data_available = "yes_series_matrix"
                    normalized_matrix_url = $read.url
                    raw_data_availability = $supp.raw_hint
                    supplementary_url = $supp.url
                    overlap_flag = $cfg.overlap_flag
                    planned_role = $cfg.role
                    story_scope = $cfg.story
                    include_primary_analysis = $include
                    exclusion_reason = $exclusion
                    curation_notes = $cfg.notes
                })
            }
            $datasetRows.Add([pscustomobject]@{
                dataset_accession = $cfg.accession
                matrix_file = $mf
                n_matrix_samples = $sampleRows.Count
                platform = $platform
                role = $cfg.role
                story_scope = $cfg.story
                include_primary_analysis = $cfg.include_primary_analysis
                normalized_matrix_url = $read.url
                supplementary_url = $supp.url
                raw_data_availability = $supp.raw_hint
                supplementary_file_count = $supp.file_count
                supplementary_examples = $supp.example_files
                overlap_flag = $cfg.overlap_flag
                notes = $cfg.notes
            })
        }
    }
    catch {
        $errors.Add([pscustomobject]@{
            dataset_accession = $cfg.accession
            error = $_.Exception.Message
            position = $_.InvocationInfo.PositionMessage
            script_stack = $_.ScriptStackTrace
        })
    }
}

$endpointMap = @{
    "GSE25066" = "pathologic_response_pcr_rd -> pCR/RD; NA excluded"
    "GSE20271" = "pcr or rd -> pCR/RD"
    "GSE32646" = "pathologic response pcr ncr -> pCR/RD where nCR maps to RD/non-pCR"
    "GSE41998" = "pcr -> Yes/No; 0 treated as missing/unmapped"
    "GSE50948" = "pcr -> pCR/RD"
    "GSE66305" = "pcr (1=yes) -> 1 pCR, 0 RD"
    "GSE163882" = "response to nac -> pCR/RD"
    "GSE194040" = "pcr -> 1 pCR, 0 RD"
    "GSE106977" = "pathological complete response -> yes pCR, no RD"
    "GSE109710" = "pcr -> 1 pCR, 0 RD"
    "GSE130786" = "drug response -> PCR/RD"
}

foreach ($cfg in $datasetConfigs) {
    $dictionaryRows.Add([pscustomobject]@{
        dataset_accession = $cfg.accession
        standard_field = "endpoint_binary_pcr1_rd0"
        source_mapping = $endpointMap[$cfg.accession]
        handling = "Unmapped or NA endpoint rows are retained in manifest but excluded from primary analysis"
    })
    $dictionaryRows.Add([pscustomobject]@{
        dataset_accession = $cfg.accession
        standard_field = "pretreatment vs on-treatment"
        source_mapping = "baseline/pretreatment inferred from series design unless a biopsy time/time point field indicates otherwise"
        handling = "Non-pretreatment samples are marked include_primary_analysis=no if detected"
    })
    $dictionaryRows.Add([pscustomobject]@{
        dataset_accession = $cfg.accession
        standard_field = "ER/PR/HER2/HR status"
        source_mapping = "dataset-specific receptor keys plus study-level defaults for TNBC/HER2-positive trials"
        handling = "Raw values preserved; normalized values are Positive/Negative/blank where possible"
    })
    $dictionaryRows.Add([pscustomobject]@{
        dataset_accession = $cfg.accession
        standard_field = "overlap_flag"
        source_mapping = $cfg.overlap_flag
        handling = "Used to prevent duplicate discovery/validation assignment"
    })
    $dictionaryRows.Add([pscustomobject]@{
        dataset_accession = $cfg.accession
        standard_field = "gsm_sample_id"
        source_mapping = "GEO !Sample_geo_accession"
        handling = "Primary sample identifier; checked for duplicate GSM IDs across the manifest"
    })
    $dictionaryRows.Add([pscustomobject]@{
        dataset_accession = $cfg.accession
        standard_field = "patient_id"
        source_mapping = "patient id/patient identifier/sample id when present; otherwise Sample_title fallback"
        handling = "Used for within-dataset repeated-patient checks; cross-GSE numeric collisions are treated as unconfirmed"
    })
    $dictionaryRows.Add([pscustomobject]@{
        dataset_accession = $cfg.accession
        standard_field = "treatment_arm"
        source_mapping = "arm/treatment arm/treatment group/preoperative treatment/type_taxane and related dataset-specific fields"
        handling = "Raw treatment-arm text is preserved; blank means not available in series matrix header"
    })
    $dictionaryRows.Add([pscustomobject]@{
        dataset_accession = $cfg.accession
        standard_field = "platform"
        source_mapping = "GPL ID inferred from matrix file name or !Sample_platform_id"
        handling = "Used to separate microarray/RNA-seq/NanoString and I-SPY2 platform-specific stress tests"
    })
    $dictionaryRows.Add([pscustomobject]@{
        dataset_accession = $cfg.accession
        standard_field = "raw_data_availability"
        source_mapping = "GEO supplementary file listing"
        handling = "Programmatic availability hint only; raw-file usability still requires dataset-specific download audit"
    })
    $dictionaryRows.Add([pscustomobject]@{
        dataset_accession = $cfg.accession
        standard_field = "include_primary_analysis"
        source_mapping = "dataset role plus endpoint/timepoint eligibility checks"
        handling = "yes = full benchmark; yes_stress_test = large platform-specific stress test; subtype_only = subtype story only; no = excluded with reason"
    })
}

$allRowsArray = @()
foreach ($r in $allRows) { $allRowsArray += $r }

$gsmDupHash = @{}
foreach ($g in @($allRowsArray | Where-Object { $_.gsm_sample_id -ne "" } | Group-Object gsm_sample_id | Where-Object { $_.Count -gt 1 })) {
    $gsmDupHash[$g.Name] = $true
}

$withinPatientDupHash = @{}
foreach ($g in @($allRowsArray | Where-Object { $_.patient_id -ne "" } | Group-Object dataset_accession,patient_id | Where-Object { $_.Count -gt 1 })) {
    $withinPatientDupHash[$g.Name] = $true
}

$crossPatientCollisionHash = @{}
foreach ($g in @($allRowsArray | Where-Object { $_.patient_id -ne "" } | Group-Object patient_id | Where-Object { ($_.Group.dataset_accession | Select-Object -Unique).Count -gt 1 })) {
    $crossPatientCollisionHash[$g.Name] = $true
}

foreach ($row in $allRowsArray) {
    $analysisKey = "$($row.dataset_accession)|$($row.matrix_file)|$($row.gsm_sample_id)"
    $withinKey = "$($row.dataset_accession), $($row.patient_id)"
    $duplicateStatus = "unique_by_gsm_and_patient_within_dataset"
    $dedupAction = "retain"

    if ($row.gsm_sample_id -ne "" -and $gsmDupHash.ContainsKey($row.gsm_sample_id)) {
        $duplicateStatus = "duplicate_gsm"
        $dedupAction = "exclude_or_resolve_before_analysis"
    }
    elseif ($row.patient_id -ne "" -and $withinPatientDupHash.ContainsKey($withinKey)) {
        $duplicateStatus = "within_dataset_patient_duplicate_or_cross_platform_repeat"
        if ($row.dataset_accession -eq "GSE194040") {
            $dedupAction = "retain_for_platform_specific_stress_tests_only; do_not_pool_I-SPY2_platforms_without_patient_level_deduplication"
        }
        else {
            $dedupAction = "review_before_patient_level_analysis"
        }
    }
    elseif ($row.patient_id -ne "" -and $crossPatientCollisionHash.ContainsKey($row.patient_id)) {
        $duplicateStatus = "cross_dataset_patient_id_collision_unconfirmed"
        $dedupAction = "retain; patient IDs are not globally unique across GSE unless proven by metadata"
    }

    $row | Add-Member -NotePropertyName analysis_unique_key -NotePropertyValue $analysisKey -Force
    $row | Add-Member -NotePropertyName duplicate_status -NotePropertyValue $duplicateStatus -Force
    $row | Add-Member -NotePropertyName dedup_action -NotePropertyValue $dedupAction -Force
}

$dedupSummary = @(
    [pscustomobject]@{ qc_item = "total_manifest_rows"; value = $allRowsArray.Count; interpretation = "sample-level rows parsed from selected GEO series matrices" },
    [pscustomobject]@{ qc_item = "duplicate_gsm_ids"; value = $gsmDupHash.Count; interpretation = "0 means no exact GSM duplicate across selected matrices" },
    [pscustomobject]@{ qc_item = "within_dataset_patient_duplicates"; value = $withinPatientDupHash.Count; interpretation = "repeated patient IDs within the same GSE; inspect before pooled patient-level analysis" },
    [pscustomobject]@{ qc_item = "cross_dataset_patient_id_collisions"; value = $crossPatientCollisionHash.Count; interpretation = "unconfirmed collisions because patient IDs are not globally unique across independent GSEs" }
)

$manifestPath = Join-Path $OutDir "analysis_ready_manifest_v1.csv"
$datasetPath = Join-Path $OutDir "dataset_level_manifest_v1.csv"
$dictPath = Join-Path $OutDir "phenotype_dictionary_v1.csv"
$errorPath = Join-Path $OutDir "manifest_build_errors_v1.csv"
$dedupPath = Join-Path $OutDir "deduplication_qc_v1.csv"
$dedupFlagPath = Join-Path $OutDir "deduplication_flagged_samples_v1.csv"

$allRows | Export-Csv -LiteralPath $manifestPath -NoTypeInformation -Encoding UTF8
$datasetRows | Export-Csv -LiteralPath $datasetPath -NoTypeInformation -Encoding UTF8
$dictionaryRows | Export-Csv -LiteralPath $dictPath -NoTypeInformation -Encoding UTF8
$dedupSummary | Export-Csv -LiteralPath $dedupPath -NoTypeInformation -Encoding UTF8
@($allRowsArray | Where-Object { $_.duplicate_status -ne "unique_by_gsm_and_patient_within_dataset" }) |
    Select-Object dataset_accession,matrix_file,gsm_sample_id,patient_id,sample_title,platform,endpoint_standard,duplicate_status,dedup_action,include_primary_analysis |
    Export-Csv -LiteralPath $dedupFlagPath -NoTypeInformation -Encoding UTF8
if ($errors.Count -gt 0) {
    $errors | Export-Csv -LiteralPath $errorPath -NoTypeInformation -Encoding UTF8
}
else {
    @([pscustomobject]@{ dataset_accession = ""; error = "" }) | Export-Csv -LiteralPath $errorPath -NoTypeInformation -Encoding UTF8
}

$summary = $allRows |
    Group-Object dataset_accession |
    ForEach-Object {
        $rows = @($_.Group)
        [pscustomobject]@{
            dataset_accession = $_.Name
            n_rows = $rows.Count
            n_endpoint_mapped = @($rows | Where-Object { $_.endpoint_binary_pcr1_rd0 -ne "" }).Count
            n_pcr = @($rows | Where-Object { $_.endpoint_binary_pcr1_rd0 -eq "1" }).Count
            n_rd = @($rows | Where-Object { $_.endpoint_binary_pcr1_rd0 -eq "0" }).Count
            include_values = ((@($rows.include_primary_analysis) | Select-Object -Unique) -join ";")
            roles = ((@($rows.planned_role) | Select-Object -Unique) -join ";")
        }
    } | Sort-Object dataset_accession

$summaryPath = Join-Path $OutDir "manifest_summary_v1.csv"
$summary | Export-Csv -LiteralPath $summaryPath -NoTypeInformation -Encoding UTF8

$summary | Format-Table -AutoSize
Write-Host "Wrote:"
Write-Host "  $manifestPath"
Write-Host "  $datasetPath"
Write-Host "  $dictPath"
Write-Host "  $summaryPath"
Write-Host "  $errorPath"
Write-Host "  $dedupPath"
Write-Host "  $dedupFlagPath"
