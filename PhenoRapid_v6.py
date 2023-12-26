# ###############################################################################################################################################################################################################################################################################
# ##################### * Load required libraries * ############################################################################################################################################################################################################################
# ###############################################################################################################################################################################################################################################################################

import os
import time
import pandas as pd
import pickle
import scipy.io
import PhenoRapid_modules_v6 as phr


# ###############################################################################################################################################################################################################################################################################
# ##################### * User Input * ##########################################################################################################################################################################################################################################
# ###############################################################################################################################################################################################################################################################################

project_path ='./Codex-Projects/Testing_Projects/'  # Edit here
project_name = 'CTLA4PD1_D10_N59_10_13_21'
# Parameters 
Resolution=1
Stardist_mask_exist=0
Avg_intensity_table_exist=0

sample_IDs=[]
sample_IDs.append('CTLA4PD1_D10_N59_10_13_21')

for sample_ID in sample_IDs:
    if not os.path.isfile(os.path.join(project_path, project_name, 'Image-Data', sample_ID, sample_ID + '_marker_info_table.csv')):
        phr.Create_Marker_Info_Table(project_path, project_name, sample_ID, includemarker='True')
    else:
        print('Marker Info Table of ' + sample_ID + ' already exists')


# Please edit /Image-Data/sample/sample_marker_info_table.csv of each sample before running next section

###############################################################################################################################################################################################################################################################################
##################### * Analysis (individual samples) * ############################################################################################################################################################################################################################################
###############################################################################################################################################################################################################################################################################
for sample_ID in sample_IDs:
    startTime = time.time()  
    print('Processing:' + sample_ID)
    qptiff_path = os.path.join('Image-Data',sample_ID)
    analysis_path = os.path.join('Analysis',sample_ID)
    CodexObj = {"project_path":[],"project_name":[],"qptiff_path":[],"sample_ID":[],"marker_info_table":[],"Labeled_Dapi_mask":[],"Labeled_Membrane_mask":[],"cell_info_table":[],"Average_Intensity_table":[],"adataScaled":[],"adata":[],"protDF":[]}
    CodexObj["project_path"]=project_path
    CodexObj["project_name"]=project_name
    CodexObj["qptiff_path"]=qptiff_path
    CodexObj["analysis_path"]=analysis_path
    CodexObj["sample_ID"]=sample_ID
    CodexObj["marker_info_table"]=pd.read_csv(os.path.join(project_path, project_name, qptiff_path, sample_ID + '_marker_info_table.csv'))
    CodexObj = phr.Stardist_segmentation(CodexObj, CLAHE=1, factor=1, tilesize=3096, min_area=5, runby='GPU', mask_already_exist=Stardist_mask_exist) #tilesize=3096, if you get an error because of min-overlap then try 2048 or 1024, at this point, it is little slower, but you will not get an error
    CodexObj = phr.Membrane_seg(CodexObj, Rin=3, Rout=15)
    CodexObj = phr.Create_cell_info_table(CodexObj)  
    CodexObj = phr.Calculate_Average_Intensity(CodexObj, AVIT=Avg_intensity_table_exist)  
    CodexObj = phr.Run_scanpy(CodexObj, runby='sample', Reso=Resolution, n_neighbors=20, min_dist=0.05, use='GPU') 
    phr.plot_CorrMatrix(CodexObj, plotby='allmarkers')  # Either 'use_of_clustering' or 'allmarkers'                  
    phr.plot_CorrMatrix(CodexObj, plotby='use_of_clustering')  # Either 'use_of_clustering' or 'allmarkers'
    phr.plot_umap(CodexObj, colorby='marker')   # Either by marker, cluster, sample, celltype
    phr.plot_umap(CodexObj, colorby='cluster', Reso=Resolution)
    phr.plot_HeatMap(CodexObj, Reso=Resolution, plotby='cluster', useonly='use_of_clustering')  # plotby: Either 'cluster' or 'celltype', useonly: Either 'use_of_clustering' or 'allmarkers'
    phr.plot_ScatterSpatial(CodexObj, colorby='cluster', Reso=Resolution)
    phr.plot_Pie_chart(CodexObj, colorby='cluster', Reso=Resolution)
    CodexObj["time_sample"]= time.time() - startTime 
    # Save Codex Object
    mypath=os.path.join(os.path.join(CodexObj["project_path"], CodexObj["project_name"], CodexObj["analysis_path"],'Output-data'))
    objectpath=mypath+'/'+'CodexObject.pkl'
    with open(objectpath, "wb") as tf:
       pickle.dump(CodexObj,tf)
    phr.Create_PPT_slides(CodexObj, runby='sample', Reso=Resolution)  # Either sample or combined
    phr.Save_CodexObj_as_Mat(CodexObj, runby='sample')
    del CodexObj   
 

# ###############################################################################################################################################################################################################################################################################
# ##################### * Analysis (combined samples) * ############################################################################################################################################################################################################################################
# ###############################################################################################################################################################################################################################################################################
# run combined dataset pipeline (if there are more than 1 sample)

if len(sample_IDs) > 1:
    startTime = time.time()  
    CODEXobj_combined = {"project_path":[],"project_name":[],"cell_info_table":[],"Average_Intensity_table":[]} 
    CODEXobj_combined["project_path"]=project_path
    CODEXobj_combined["project_name"]=project_name
    CODEXobj_combined["analysis_path"]=os.path.join('Analysis','Combined_dataset')
    CODEXobj_combined["sample_IDs"]=sample_IDs
    print('Start analyzing all the dataset samples combined :)')
    CODEXobj_combined = phr.Run_dataset_combined(CODEXobj_combined)
    CODEXobj_combined = phr.Run_scanpy(CODEXobj_combined, runby='combined' ,Reso=Resolution, n_neighbors=20, min_dist=0.05)
    phr.plot_umap(CODEXobj_combined, colorby='sample')
    phr.plot_Pie_chart(CODEXobj_combined, colorby='sample')
    phr.plot_CorrMatrix(CODEXobj_combined, plotby='allmarkers')  # Either 'use_of_clustering' or 'allmarkers'    # I put allmarkers here because Run_dataset_combined has already choose only the use for clustering
    phr.plot_umap(CODEXobj_combined, colorby='marker') 
    phr.plot_umap(CODEXobj_combined, colorby='cluster', Reso=Resolution) 
    phr.plot_HeatMap(CODEXobj_combined, Reso=Resolution, plotby='cluster', useonly='allmarkers')
    phr.plot_Pie_chart(CODEXobj_combined, colorby='cluster', Reso=Resolution)
    phr.plot_ScatterSpatial(CODEXobj_combined, colorby='sample', Reso=Resolution)
    phr.Create_PPT_slides(CODEXobj_combined, runby='combined', Reso=Resolution)
    CODEXobj_combined["time_Combined"]= time.time() - startTime 
    phr.Save_CodexObj_as_Mat(CODEXobj_combined, runby='combined')   
    del CODEXobj_combined
    



# ###############################################################################################################################################################################################################################################################################
# ##################### * Run pipeline after phenotyping * ############################################################################################################################################################################################################################################
# ###############################################################################################################################################################################################################################################################################

# for sample_ID in sample_IDs:
    
#     startTime = time.time()  
#     print('Processing:' + sample_ID)
#     analysis_path = os.path.join('Analysis',sample_ID)
#     mypath=os.path.join(os.path.join(project_path, project_name, analysis_path,'Output-data'))
#     objectpath=mypath+'/'+'CodexObject.pkl'
#     with open(objectpath, "rb") as tf:
#         CodexObj = pickle.load(tf)
    
#     phr.prepare_pheno_objects(CodexObj, Reso=Resolution)
#     phr.plot_umap(CodexObj, colorby='celltype', Reso=Resolution)
#     phr.plot_HeatMap(CodexObj,Reso=Resolution, plotby='celltype', useonly='use_of_clustering')  # Either 'use_of_clustering' or 'celltype'
#     phr.plot_ScatterSpatial(CodexObj, colorby='celltype', Reso=Resolution)
#     phr.plot_Pie_chart(CodexObj, colorby='celltype', Reso=Resolution)
#     phr.Create_PPT_slides(CodexObj, runby='celltype', Reso=Resolution)
#     CodexObj["time_sample"]= time.time() - startTime 
#     # phr.Save_CodexObj_as_Mat(CodexObj, runby='sample')
#     del CodexObj   



# objectpath=mypath+'/'+'CodexObject.pkl'
# with open(objectpath, "rb") as tf:
#    pickle.load(tf)  





