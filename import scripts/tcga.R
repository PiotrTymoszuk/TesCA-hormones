# Import of the TCGA study data, version by Liu 2018 deposited at cBioportal

  insert_head()

# container -------

  tcga <- list()

# clinical information -------

  insert_msg('Clinical information')

  ## patient

  tcga$clinic$patient <- read_tsv('./data/TCGA/data_clinical_patient.txt',
                                  skip = 4)

  tcga$clinic$patient <- tcga$clinic$patient %>%
    transmute(patient_id = PATIENT_ID,
              histology = stri_replace(SUBTYPE,
                                       regex = '^.*_', replacement = ''),
              histology = car::recode(histology, "'non-seminoma' = 'NSGCT'"),
              histology = factor(histology, c('seminoma', 'NSGCT')),
              histology_icd = car::recode(ICD_O_3_HISTOLOGY,
                                             "'9061/3' = 'seminoma';
                                             '9070/3' = 'embryonal carcinoma';
                                             '9071/3' = 'yolk sac cancer';
                                             '9080/0' = 'benign teratoma';
                                             '9080/3' = 'malignant teratoma';
                                             '9081/3' = 'teratocarcinoma';
                                             '9085/3' = 'germinal mixed histology'"),
              histology_icd = factor(histology_icd),
              histology_icd = fct_relevel(histology_icd,
                                          'seminoma',
                                          'germinal mixed histology'),
              age = as.numeric(AGE),
              pt_stage = stri_replace(AJCC_PATHOLOGIC_TUMOR_STAGE,
                                      regex = '.*\\s{1}',
                                      replacement = ''),
              pt_stage = stri_replace(pt_stage,
                                      regex = '(A|B|S)$',
                                      replacement = ''),
              pt_stage = factor(pt_stage, c('I', 'II', 'III', 'IV')),
              fup_days = as.numeric(DAYS_LAST_FOLLOWUP),
              neoadjuvant = tolower(HISTORY_NEOADJUVANT_TRTYN),
              neoadjuvant = factor(neoadjuvant, c('no', 'yes')),
              pm_stage = stri_extract(PATH_M_STAGE, regex = '^M\\d{1}'),
              pm_stage = factor(pm_stage, c('M0', 'M1')),
              pn_stage = stri_extract(PATH_N_STAGE,
                                      regex = '^N\\d{1}'),
              pn_stage = factor(pn_stage, c('N0', 'N1', 'N2')),
              path_t_stage = stri_extract(PATH_T_STAGE,
                                          regex = '^T\\d{1}'),
              path_t_stage = factor(path_t_stage, c('T1', 'T2', 'T3', 'T4')),
              race = factor(RACE),
              radiotherapy = tolower(RADIATION_THERAPY),
              radiotherapy = factor(radiotherapy, c('no', 'yes')),
              death = stri_extract(OS_STATUS, regex = '^\\d{1}'),
              death = as.numeric(death),
              os_days = as.numeric(OS_MONTHS) * 30.437,
              tumor_death = stri_extract(DSS_STATUS, regex = '^\\d{1}'),
              tumor_death = as.numeric(tumor_death),
              dss_days = as.numeric(DSS_MONTHS) * 30.437,
              relapse = stri_extract(DFS_STATUS, regex = '^\\d{1}'),
              relapse = as.numeric(relapse),
              rfs_days = as.numeric(DFS_MONTHS) * 30.437,
              progression = stri_extract(PFS_STATUS, regex = '^\\d{1}'),
              progression = as.numeric(progression),
              progression_factor = ifelse(progression == 1, 'yes', 'no'),
              progression_factor = factor(progression_factor, c('no', 'yes')),
              pfs_days = as.numeric(PFS_MONTHS) * 30.437,
              relapse = ifelse(is.na(relapse) & !is.na(progression_factor),
                               ifelse(progression_factor == 'yes',
                                      1, 0),
                               relapse),
              rfs_days = ifelse(is.na(rfs_days), pfs_days, rfs_days),
              relapse_factor = ifelse(relapse == 1, 'yes', 'no'),
              relapse_factor = factor(relapse_factor, c('no', 'yes')),
              death_factor = ifelse(death == 1, 'yes', 'no'),
              death_factor = factor(death_factor, c('no', 'yes')))

  ## sample

  tcga$clinic$sample <- read_tsv('./data/TCGA/data_clinical_sample.txt',
                                  skip = 4)

  tcga$clinic$sample <- tcga$clinic$sample %>%
    transmute(sample_id = SAMPLE_ID,
              patient_id = PATIENT_ID,
              oncotree = factor(ONCOTREE_CODE),
              oncotree = fct_relevel(oncotree, 'SEM', 'MGCT', 'EMBCA'),
              histology_detail = factor(TUMOR_TYPE),
              aneuploidy_score = as.numeric(ANEUPLOIDY_SCORE),
              mantis_msi_score = as.numeric(MSI_SCORE_MANTIS),
              sensor_msi_score = as.numeric(MSI_SENSOR_SCORE),
              tmb = as.numeric(TMB_NONSYNONYMOUS))

  ## risk stratification at diagnosis

  tcga$clinic$risk <-
    read_tsv('./data/TCGA/data_timeline_status.txt') %>%
    filter(STATUS == 'Initial Diagnosis') %>%
    transmute(patient_id = PATIENT_ID,
              marker_status = toupper(SERUM_MARKERS),
              marker_status = factor(marker_status,
                                     c('S0', 'S1', 'S2', 'S3')),
              IGCCCG_risk_group = factor(tolower(IGCCCG_STAGE),
                                         c('good', 'intermediate', 'poor')))

  ## merging the clinical data together

  tcga$clinic <- tcga$clinic %>%
    reduce(left_join, by = 'patient_id') %>%
    relocate(sample_id, patient_id)

# New tumor events ---------

  insert_msg('New tumor events')

  ## appending with the sample ID as well

  tcga$new_event <-
    read_tsv('./data/TCGA/data_timeline_status.txt') %>%
    filter(STATUS %in% c('Locoregional Recurrence',
                         'New Primary Tumor',
                         'Distant Metastasis',
                         'Last Follow Up')) %>%
    transmute(patient_id = PATIENT_ID,
              event_type = STATUS,
              location = ANATOMIC_SITE,
              event_days = as.numeric(START_DATE)) %>%
    left_join(tcga$clinic[c('patient_id', 'sample_id')],
              by = 'patient_id')

# Expression and expression data annotation -------

  insert_msg('Expression and expression annotation')

  ## keeps genes with symbols. Manual correction/check
  ## of the duplicated symbols

  tcga$expression <- read_tsv('./data/TCGA/data_mrna_seq_v2_rsem.txt') %>%
    mutate(Entrez_Gene_Id = as.character(Entrez_Gene_Id),
           probe_id = paste0('probe_', 1:nrow(.)))

  tcga$annotation <- tcga$expression %>%
    transmute(probe_id = probe_id,
              entrez_id = Entrez_Gene_Id) %>%
    mutate(gene_symbol = mapIds(org.Hs.eg.db,
                                keys = entrez_id,
                                keytype = 'ENTREZID',
                                column = 'SYMBOL')) %>%
    filter(complete.cases(.))

  ## integration of expression: duplicated entries are collapsed
  ## by arithmetic mean. Log2(x + 1) transformation

  tcga$expression <- tcga$expression %>%
    select(-Entrez_Gene_Id, -Hugo_Symbol) %>%
    column_to_rownames('probe_id') %>%
    integrate_expression(annotation = tcga$annotation,
                         trans_fun = function(x) log2(x + 1))

  tcga$annotation <- tcga$annotation %>%
    filter(!duplicated(gene_symbol))

# Protein expression ----------

  insert_msg('Protein expression and its annotation')

  ## the protein lexicon is appended with proein classification
  ## made pe hand in a separate Excel file

  tcga$protein <- read_tsv('./data/TCGA/data_rppa.txt')

  tcga$protein_annotation$lexicon <-
    tibble(variable = tcga$protein$Composite.Element.REF) %>%
    mutate(protein_symbols = stri_split_fixed(variable,
                                              pattern = '|',
                                              simplify = TRUE)[, 1])

  tcga$protein_annotation$classification <-
    read_xlsx('./data/prot_lexicon.xlsx')

  tcga$protein_annotation <- reduce(tcga$protein_annotation,
                                    left_join, by = 'variable')

  ## the protein expression levels are already at the log2-scale

  tcga$protein <- tcga$protein %>%
    column_to_rownames('Composite.Element.REF') %>%
    t %>%
    as.data.frame %>%
    rownames_to_column('sample_id') %>%
    as_tibble

# Mutation data --------

  insert_msg('Somatic mutation data')

  ## mutations with consequences for protein sequences

  tcga$mutation_details <- read_tsv('./data/TCGA/data_mutations.txt') %>%
    filter(Variant_Classification != 'Silent',
           !is.na(Protein_position))

  ## 0/1 coded absence/presence of mutations
  ## appending with samples without any mutations

  tcga$mutations <- tcga$mutation_details %>%
    extract_mutations

  tcga$mutations <-
    full_rbind(tcga$mutations,
               filter(tcga$clinic,
                      !sample_id %in% tcga$mutations$sample_id)['sample_id']) %>%
    map_dfc(~ifelse(is.na(.x), 0, .x))

# Gene amplifications and deletions --------

  insert_msg('Gene amplifications and deletions')

  tcga$cna <- read_tsv('./data/TCGA/data_cna.txt') %>%
    filter(!duplicated(Hugo_Symbol)) %>%
    select(-Entrez_Gene_Id, -Cytoband) %>%
    column_to_rownames('Hugo_Symbol') %>%
    t

  ## homozygous deletions and high-level amplification
  ## are included in the analyses

  tcga$deletions <- ifelse(tcga$cna < -1, 1, 0)
  tcga$amplifications <- ifelse(tcga$cna > 1, 1, 0)

  tcga[c("deletions", "amplifications")] <-
    tcga[c("deletions", "amplifications")] %>%
    map(as.data.frame) %>%
    map(rownames_to_column, 'sample_id') %>%
    map(as_tibble)

# Therapy data -------

  insert_msg('Therapy information')

  tcga$treatment <- read_tsv('./data/TCGA/data_timeline_treatment.txt')

  tcga$treatment <- tcga$treatment %>%
    transmute(patient_id = PATIENT_ID,
              treatment_type = factor(tolower(TREATMENT_TYPE),
                                      c('chemotherapy',
                                        'hormone therapy',
                                        'radiation therapy')),
              agent = factor(tolower(AGENT)),
              treatment_class = fct_collapse(agent,
                                             bleomycin = 'bleomycin',
                                             platinum = c('carboplatin',
                                                          'cisplatin',
                                                          'platinum',
                                                          'gemcitabine + oxaliplatin'),
                                             etoposide = 'etoposide',
                                             radiation = 'radiation 1'),
              treatment_class = fct_lump_min(treatment_class, 5,
                                             other_level = 'other'),
              response = factor(tolower(MEASURE_OF_RESPONSE)),
              response = fct_collapse(response,
                                      CR = 'complete response',
                                      PR = 'partial response',
                                      `SD/PD` = c('stable disease',
                                                  'clinical progressive disease')),
              response = fct_relevel(response,
                                     'CR', 'PR', 'SD/PD')) %>%
    left_join(tcga$clinic[, c("patient_id", "sample_id")],
              by = 'patient_id') %>%
    relocate(sample_id, patient_id)

# Caching the results ---------

  insert_msg('Caching the cleared data sets')

  tcga <- tcga[c("clinic", "new_event",
                 "expression", "annotation",
                 "protein", "protein_annotation",
                 "mutations", "deletions", "amplifications",
                 "treatment")]

  ## adding tissue information

  tcga[c("clinic",
         "expression",
         "protein",
         "mutations",
         "deletions",
         "amplifications",
         "treatment")] <-
    tcga[c("clinic",
           "expression",
           "protein",
           "mutations",
           "deletions",
           "amplifications",
           "treatment")] %>%
    map(mutate,
        tissue = factor('tumor', c('normal', 'tumor')))

  save(tcga, file = './data/tcga.RData')

# END -------

  insert_tail()
