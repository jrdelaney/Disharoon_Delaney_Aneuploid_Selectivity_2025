import os
import pandas as pd
import requests
import json
import csv

final_files_dir = 'final_TCGA_datasets'
os.makedirs(final_files_dir, exist_ok=True)

### DOWNLOAD AND COMBINE DATA
## TCGA

# Automate pulling these for all studies

data_files = [
    '1p_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '1q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '2p_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '2q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '3p_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '3q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '4p_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '4q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '5p_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '5q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '6p_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '6q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '7p_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '7q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '8p_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '8q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '9p_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '9q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '10p_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '10q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '11p_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '11q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '12p_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '12q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '13q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '14q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '15q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '16p_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '16q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '17p_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '17q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '18p_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '18q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '19p_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '19q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '20p_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '20q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '21q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
    '22q_status__Putative_arm-level_copy-number_from_GISTIC.full.txt',
]

# Initialize an empty DataFrame for combined aneuploidy data
combined_aneuploidy_df = pd.DataFrame(columns=['Study ID', 'Patient ID', 'Sample ID'])

# Iterate through each file and concatenate the data
count = 0
for file in data_files:
    temp_aneuploidy_df = pd.read_csv(os.path.join('GISTIC_pan_cancer', file), sep='\t')
    count += 1
    if count == 1:
        combined_aneuploidy_df = temp_aneuploidy_df
    else:
        combined_aneuploidy_df = combined_aneuploidy_df.merge(temp_aneuploidy_df, how='left', on=['Study ID', 'Patient ID', 'Sample ID'])

# Save the combined aneuploidy data to CSV
combined_aneuploidy_df.to_csv(os.path.join(final_files_dir, '1_cbioportal_GISTIC_aneuploidy_data_combined.csv'), index=False)
print(combined_aneuploidy_df)

# Read cancer type data
cancer_type_df = pd.read_csv(os.path.join('GISTIC_pan_cancer', 'Cancer_Type_chrm_arm_GISTIC.txt'), sep='\t')

# Merge aneuploidy data with cancer type information
merged_aneuploidy_cancer_df = combined_aneuploidy_df.merge(cancer_type_df, how='left', on='Sample ID')

# Save the merged DataFrame
merged_aneuploidy_cancer_df.to_csv(os.path.join(final_files_dir, '2_cbioportal_GISTIC_aneuploidy_data_combined_with_cancer_type.csv'), index=False)

# Rename long columns
rename_dict = {
    '1p_status: Putative arm-level copy-number from GISTIC': '1p',
    '1q_status: Putative arm-level copy-number from GISTIC': '1q',
    '2p_status: Putative arm-level copy-number from GISTIC': '2p',
    '2q_status: Putative arm-level copy-number from GISTIC': '2q',
    '3p_status: Putative arm-level copy-number from GISTIC': '3p',
    '3q_status: Putative arm-level copy-number from GISTIC': '3q',
    '4p_status: Putative arm-level copy-number from GISTIC': '4p',
    '4q_status: Putative arm-level copy-number from GISTIC': '4q',
    '5p_status: Putative arm-level copy-number from GISTIC': '5p',
    '5q_status: Putative arm-level copy-number from GISTIC': '5q',
    '6p_status: Putative arm-level copy-number from GISTIC': '6p',
    '6q_status: Putative arm-level copy-number from GISTIC': '6q',
    '7p_status: Putative arm-level copy-number from GISTIC': '7p',
    '7q_status: Putative arm-level copy-number from GISTIC': '7q',
    '8p_status: Putative arm-level copy-number from GISTIC': '8p',
    '8q_status: Putative arm-level copy-number from GISTIC': '8q',
    '9p_status: Putative arm-level copy-number from GISTIC': '9p',
    '9q_status: Putative arm-level copy-number from GISTIC': '9q',
    '10p_status: Putative arm-level copy-number from GISTIC': '10p',
    '10q_status: Putative arm-level copy-number from GISTIC': '10q',
    '11p_status: Putative arm-level copy-number from GISTIC': '11p',
    '11q_status: Putative arm-level copy-number from GISTIC': '11q',
    '12p_status: Putative arm-level copy-number from GISTIC': '12p',
    '12q_status: Putative arm-level copy-number from GISTIC': '12q',
    '13q_status: Putative arm-level copy-number from GISTIC': '13q',
    '14q_status: Putative arm-level copy-number from GISTIC': '14q',
    '15q_status: Putative arm-level copy-number from GISTIC': '15q',
    '16p_status: Putative arm-level copy-number from GISTIC': '16p',
    '16q_status: Putative arm-level copy-number from GISTIC': '16q',
    '17p_status: Putative arm-level copy-number from GISTIC': '17p',
    '17q_status: Putative arm-level copy-number from GISTIC': '17q',
    '18p_status: Putative arm-level copy-number from GISTIC': '18p',
    '18q_status: Putative arm-level copy-number from GISTIC': '18q',
    '19p_status: Putative arm-level copy-number from GISTIC': '19p',
    '19q_status: Putative arm-level copy-number from GISTIC': '19q',
    '20p_status: Putative arm-level copy-number from GISTIC': '20p',
    '20q_status: Putative arm-level copy-number from GISTIC': '20q',
    '21q_status: Putative arm-level copy-number from GISTIC': '21q',
    '22q_status: Putative arm-level copy-number from GISTIC': '22q',
}

# Rename the chromosome arm columns
merged_aneuploidy_cancer_df = merged_aneuploidy_cancer_df.rename(columns=rename_dict)

# Save the renamed DataFrame
merged_aneuploidy_cancer_df.to_csv(os.path.join(final_files_dir, '3_cbioportal_GISTIC_aneuploidy_data_combined_with_cancer_type.csv'), index=False)

# Read combined clinical data
clinical_data_df = pd.read_csv(os.path.join('GISTIC_pan_cancer', 'combined_study_clinical_data.tsv'), sep='\t')

# Merge with the aneuploidy and cancer type data
merged_full_df = merged_aneuploidy_cancer_df.merge(clinical_data_df, how='left', left_on =['Study ID_x','Patient ID_x','Sample ID'] ,right_on=['Study ID', 'Patient ID', 'Sample ID'])

# Save the fully merged DataFrame
merged_full_df.to_csv(os.path.join(final_files_dir, '4_cbio_aneuploidy_data.csv'), index=False)

# Setup the endpoint for making cbioportal data requests
all_study_link = 'https://www.cbioportal.org/study/summary?id=paac_jhu_2014%2Cmel_tsam_liang_2017%2Call_stjude_2016%2Caml_ohsu_2022%2Caml_ohsu_2018%2Claml_tcga_pan_can_atlas_2018%2Cmnm_washu_2016%2Cacyc_fmi_2014%2Cacyc_jhu_2016%2Cacyc_mda_2015%2Cacyc_mskcc_2013%2Cacyc_sanger_2013%2Cacc_2019%2Cacbc_mskcc_2015%2Cacc_tcga_pan_can_atlas_2018%2Campca_bcm_2016%2Cbcc_unige_2016%2Cbladder_columbia_msk_2018%2Cblca_mskcc_solit_2014%2Cblca_mskcc_solit_2012%2Cblca_bgi%2Cblca_dfarber_mskcc_2014%2Cblca_tcga_pan_can_atlas_2018%2Clgg_tcga_pan_can_atlas_2018%2Cbrca_hta9_htan_2022%2Cbrca_metabric%2Cbrca_mskcc_2019%2Cbrca_smc_2018%2Cbfn_duke_nus_2015%2Cbrca_bccrc%2Cbrca_broad%2Cbrca_sanger%2Cbrca_tcga_pan_can_atlas_2018%2Ccesc_tcga_pan_can_atlas_2018%2Cpan_origimed_2020%2Cchol_icgc_2017%2Cchol_nccs_2013%2Cchol_nus_2012%2Cchol_tcga_pan_can_atlas_2018%2Clcll_broad_2013%2Ccll_broad_2015%2Ccll_broad_2022%2Ccll_iuopa_2015%2Ccllsll_icgc_2011%2Cccrcc_dfci_2019%2Ccoad_caseccc_2015%2Ccoad_cptac_2019%2Ccoad_silu_2022%2Ccoadread_dfci_2016%2Ccoadread_genentech%2Ccoadread_tcga_pan_can_atlas_2018%2Ccoadread_mskcc%2Ccoadread_cass_2020%2Chccihch_pku_2019%2Ccscc_dfarber_2015%2Ccscc_hgsc_bcm_2014%2Ccscc_ucsf_2021%2Cctcl_columbia_2015%2Cpact_jhu_2011%2Cdesm_broad_2015%2Cdifg_glass_2019%2Cdlbcl_dfci_2018%2Cdlbcl_duke_2017%2Cdlbc_tcga_pan_can_atlas_2018%2Cnhl_bcgsc_2013%2Ccrc_nigerian_2020%2Cucec_cptac_2020%2Cucec_ccr_msk_2022%2Cucec_ccr_cfdna_msk_2022%2Cesca_broad%2Cesca_tcga_pan_can_atlas_2018%2Cescc_icgc%2Cescc_ucla_2014%2Ces_iocurie_2014%2Cgbc_shanghai_2014%2Cegc_tmucih_2015%2Cstad_oncosg_2018%2Cgbm_cptac_2021%2Cgbm_columbia_2019%2Cgbm_tcga_pan_can_atlas_2018%2Cglioma_msk_2018%2Chnsc_broad%2Chnsc_jhu%2Chnsc_tcga_pan_can_atlas_2018%2Cliad_inserm_fr_2014%2Chcc_meric_2021%2Chcc_inserm_fr_2015%2Chistiocytosis_cobi_msk_2019%2Cpanet_shanghai_2013%2Cchol_jhu_2013%2Cihch_ismms_2015%2Cihch_smmu_2014%2Ckich_tcga_pan_can_atlas_2018%2Ckirc_bgi%2Cccrcc_irc_2014%2Ckirc_tcga_pan_can_atlas_2018%2Ckirp_tcga_pan_can_atlas_2018%2Chcc_msk_venturaa_2018%2Clihc_amc_prv%2Clihc_riken%2Clihc_tcga_pan_can_atlas_2018%2Clgg_ucsf_2014%2Clgsoc_mapk_msk_2022%2Cluad_broad%2Cluad_cptac_2020%2Cluad_oncosg_2020%2Cluad_tcga_pan_can_atlas_2018%2Cluad_tsp%2Clung_smc_2016%2Clung_nci_2022%2Clusc_cptac_2021%2Clusc_tcga_pan_can_atlas_2018%2Cmsk_impact_2017%2Cmixed_allen_2018%2Cmpnst_mskcc%2Cmcl_idibips_2013%2Cmbl_broad_2012%2Cmbl_dkfz_2017%2Cmbl_pcgp%2Cmbl_sickkids_2016%2Cskcm_mskcc_2014%2Cmng_utoronto_2021%2Cmeso_tcga_pan_can_atlas_2018%2Cbiliary_tract_summit_2022%2Cbrca_igr_2015%2Cmel_dfci_2019%2Cskcm_dfci_2015%2Cskcm_vanderbilt_mskcc_2015%2Cmel_ucla_2016%2Cprad_mich%2Cprad_su2c_2019%2Cmetastatic_solid_tumors_mich_2017%2Cmixed_selpercatinib_2020%2Cmm_broad%2Cmds_tokyo_2011%2Cmds_iwg_2022%2Cmpn_cimr_2013%2Cstmyec_wcm_2022%2Cnpc_nusingapore%2Cnbl_amc_2012%2Cnbl_ucologne_2015%2Cnepc_wcm_2016%2Cnhl_bcgsc_2011%2Cnsclc_mskcc_2018%2Cnsclc_tracerx_2017%2Cnsclc_unito_2016%2Chnsc_mdanderson_2013%2Cov_tcga_pan_can_atlas_2018%2Cpog570_bcgsc_2020%2Cpancan_pcawg_2020%2Cpaad_qcmg_uq_2016%2Cpaad_tcga_pan_can_atlas_2018%2Cpaad_utsw_2015%2Cpaad_cptac_2021%2Cpanet_jhu_2011%2Cpanet_arcnet_2017%2Call_phase2_target_2018_pub%2Caml_target_2018_pub%2Cbrain_cptac_2020%2Ces_dfarber_broad_2014%2Cnbl_target_2018_pub%2Cpediatric_dkfz_2017%2Cmixed_pipseq_2017%2Cpptc_2019%2Crt_target_2018_pub%2Cwt_target_2018_pub%2Cpcpg_tcga_pan_can_atlas_2018%2Cplmeso_nyu_2015%2Ccrc_hta11_htan_2021%2Cpcnsl_mayo_2015%2Cprad_broad%2Cprad_fhcrc%2Cprad_mskcc%2Cprad_eururol_2017%2Cprad_tcga_pan_can_atlas_2018%2Cprad_mskcc_cheny1_organoids_2014%2Cprostate_dkfz_2018%2Cprad_msk_2019%2Cprostate_pcbm_swiss_2019%2Cprad_msk_mdanderson_2023%2Cbrca_cptac_2020%2Cccrcc_utokyo_2013%2Cnccrcc_genentech_2014%2Cmrt_bcgsc_2016%2Crms_nih_2014%2Csummit_2018%2Csarc_mskcc%2Csarc_tcga_pan_can_atlas_2018%2Cskcm_broad%2Cskcm_tcga_pan_can_atlas_2018%2Cskcm_yale%2Cskcm_broad_brafresist_2012%2Cscco_mskcc%2Csclc_jhu%2Csclc_ucologne_2015%2Csclc_cancercell_gardner_2017%2Csarcoma_msk_2022%2Cvsc_cuk_2018%2Cstad_pfizer_uhongkong%2Cstad_tcga_pan_can_atlas_2018%2Cstad_utokyo%2Ctgct_tcga_pan_can_atlas_2018%2Cangs_painter_2020%2Cangs_project_painter_2018%2Cbrca_mbcproject_wagle_2017%2Cmpcproject_broad_2021%2Ctet_nci_2014%2Cthym_tcga_pan_can_atlas_2018%2Cthca_tcga_pan_can_atlas_2018%2Curcc_mskcc_2016%2Cutuc_mskcc_2015%2Cutuc_cornell_baylor_mdacc_2019%2Cutuc_igbmc_2021%2Cutuc_msk_2019%2Cblca_bcan_hcrn_2022%2Cblca_cornell_2016%2Cucs_jhu_2014%2Cucs_tcga_pan_can_atlas_2018%2Cuccc_nih_2017%2Cucec_tcga_pan_can_atlas_2018%2Cum_qimr_2016%2Cuvm_tcga_pan_can_atlas_2018'
study1 = all_study_link.replace('https://www.cbioportal.org/study/summary?id=', "")
study_list = all_study_link.split(sep='%2')

### Get data
url = 'https://www.cbioportal.org/api/treatments/sample'
headers = {
    'Accept': 'application/json',
    'Content-Type': 'application/json',
    'Origin': 'https://www.cbioportal.org',
}

data = {
    "studyIds": study_list,
    "alterationFilter": {
        "copyNumberAlterationEventTypes": {"AMP": True, "HOMDEL": True},
        "mutationEventTypes": {"any": True},
        "structuralVariants": None,
        "includeDriver": True,
        "includeVUS": True,
        "includeUnknownOncogenicity": True,
        "includeUnknownTier": True,
        "includeGermline": True,
        "includeSomatic": True,
        "includeUnknownStatus": True,
        "tiersBooleanMap": {}
    }
}

response = requests.post(url, headers=headers, json=data)
with open('treatment.json', 'w') as json_file:
    json.dump(response.json(), json_file, indent=4)

# Open CSV file to write data
with open('treatment_data.csv', mode='w', newline='') as file:
    writer = csv.writer(file)

    # Write the header
    writer.writerow(['treatment', 'patientId', 'sampleId', 'studyId'])

    # Iterate over the JSON data and extract required fields
    for item in response.json():
        treatment = item['treatment']
        for sample in item['samples']:
            patientId = sample['patientId']
            sampleId = sample['sampleId']
            studyId = sample['studyId']

            # Write the row to CSV
            writer.writerow([treatment, patientId, sampleId, studyId])

# Combine CSV with main data set joining on sample_id
cbio_aneuploidy_df = merged_full_df
treatment_data_df = pd.read_csv('treatment_data.csv')

# Correct treatment naming
treatment_data_df['treatment'] = treatment_data_df['treatment'].str.lower()
# Split multiple treatments into lists
treatment_data_df['treatment'] = treatment_data_df['treatment'].str.replace(' + ', ', ')
treatment_data_df['treatment'] = treatment_data_df['treatment'].str.split(', ')

# Merge aneuploidy data with treatment data
merged_with_treatment_df = cbio_aneuploidy_df.merge(
    treatment_data_df,
    left_on=['Study ID', 'Patient ID', 'Sample ID'],
    right_on=['studyId','patientId','sampleId'],
    how='left'
)

# Save the merged DataFrame with treatment information
merged_with_treatment_df.to_csv(os.path.join(final_files_dir, '5_full_cbio_data_treatment.csv'), index=False)

# Drop unnecessary columns
merged_with_treatment_df = merged_with_treatment_df.drop(columns=['Cancer Type_y',])

# Rename columns for clarity
merged_with_treatment_df.rename(columns={'Cancer Type_x': 'Cancer Type'}, inplace=True)

merged_with_treatment_df.to_csv(os.path.join(final_files_dir, '6_full_clean_cbio_data_treatment.csv'), index=False)

# explode the dataframe on treatment lists
merged_with_treatment_df = merged_with_treatment_df.explode('treatment').reset_index(drop=True)

merged_with_treatment_df['treatment'] = merged_with_treatment_df['treatment'].str.replace(r'radiation\s*\d+', 'radiation', regex=True)

collapsed_df = merged_with_treatment_df.drop_duplicates()

# Now, you can save the DataFrame without altering the aneuploidy statuses
collapsed_df.to_csv(os.path.join(final_files_dir, '7_TCGA_final_exploded_treatment_aneuploidy_by_treatment.csv'), index=False)

print('Done!')