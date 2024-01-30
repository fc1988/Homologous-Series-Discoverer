while True:

    import numpy as np; import pandas as pd;
    import os; import copy
    from pyteomics import mass
    import PySimpleGUI as sg

    # %% GUI & Input
    sg.theme('Reddit')

    sample_col = [[sg.Text('Filepath'), sg.In(size=(35,1), enable_events=True ,key='-FOLDER_SAM-',tooltip = 'Input .csv file with two variables including mass and rt'), sg.FileBrowse()]]
    # create layout
    layout1 = [[sg.Text('Homologous Series Discoverer', size=(24, 1), font=('Arial',22), text_color = 'Black')],
              [sg.Text("Sample filename (.csv)")],[sg.Column(sample_col, element_justification='c')],
              [sg.Text('Homologous series differences (Space separated)')], [sg.InputText('CF2',key='frags', tooltip ='Input fragment mass differences here, either as chemical formula or accurate mass (e.g CF2). \nExample: (CF2)n, (C2F4)n, HF, CnH3F2n-3, CnH2F2n-4, CnHF2n-5, C2F4-HF = 79.98739, C2F4-H2O = 81.98304, CnF2nSO3, CF2O, CF3')],
              [sg.Text("Mass tolerance for homologous series (Da)")], [sg.InputText('0.0025', key='hs_tol')],
              [sg.Text("Minimum number of homologues")], [sg.InputText('3', key='n_min')]]

    tabgrp = [[sg.TabGroup([[sg.Tab('Mass difference', layout1)]])],[sg.Submit(tooltip='HSD\u0394S')]]

    full_path = False
    # Create the GUI window
    window = sg.Window("Homologous Series Discoverer v1.0", tabgrp, resizable = True, finalize = True)
    while True:
         [event, values] = window.read()  
         if event in (None, 'Exit'):
             break
         elif event == "Submit":
             window.close()
             # =====================================
             saved_values = values
             full_path =  values['-FOLDER_SAM-']
               
             if not values.get('frags').split():
                frags = ['CF2']
             else:
                frags = values.get('frags').split()
                for n in range(len(frags)):
                    if frags[n][0].isnumeric():
                        frags[n] = float(frags[n])
                        
             if not values.get('hs_tol'):  
                 hs_tol = 0.0025
             else:
                 hs_tol = float(values.get('hs_tol'))
                 
             if not values.get('n_min'):  
                 n_min = 3
             else:
                 n_min = float(values.get('n_min'))

    # %% Folder generation
    # =========================================================================
    # check if path is specified, otherwise end program (also triggered by close button)
    if full_path == False:
        break

    print('\nStarted data evaluation...\n...')

    # =========================================================================
    # create results folder with sample name attached using os packages
    sample_name = os.path.splitext(os.path.basename(full_path))[0]

    Res_folder = 'Results_HSD_' + sample_name

    # check if folder exists, if yes overwrite it
    if not os.path.exists(Res_folder):
        os.makedirs(Res_folder)

    # =========================================================================
    # create mass of difference array from chemical formulae

    frags = ['CF2']
    m_frags = np.array(np.zeros((len(frags))))
    for n in range(len(frags)):
      if isinstance(frags[n], str):
        m_frags[n] = mass.calculate_mass(formula = frags[n])
      else:
        m_frags[n] = frags[n]

    # convert fragment masses to string for column names
    frags_str = [[None] for x in range(len(frags))]
    for n in range(len(frags)):
      frags_str[n] = str(frags[n])

    # =========================================================================
    # read mass and rt on list
    [_, extension] = os.path.splitext(full_path)
    data = pd.read_csv(full_path)
    mz_vec_corr= data['mass']
    RT_vec_corr= data['rt']
    mz_vec_corr = np.array(mz_vec_corr)
    RT_vec_corr = np.array(RT_vec_corr)

    # =========================================================================
    # Create homologous series (HS) based on the user input
    # =========================================================================
    mz_HS_tot = np.array([])

    for j in range(len(m_frags)):
      rep_unit = m_frags[j]

    #Calculate modulo: Compounds from the same homologous series bear an identical modulo
    modulo = mz_vec_corr % rep_unit

    # Calculation of Kendrick masses and features that are in the same homologous series
    KM = mz_vec_corr*round(rep_unit)/rep_unit
    KM_round = np.round_(KM)
    KMD = KM - KM_round
    Mod_HS = {'mz':mz_vec_corr,'RT':RT_vec_corr,'mod':modulo,'KMD':KMD}
    Mod_HS_Dataframe = pd.DataFrame(data=Mod_HS)
    Mod_HS_Dataframe.sort_values('mod', inplace=True, ignore_index=True)
    HS_num = np.zeros(len(mz_vec_corr))
    exit_loop = 0

    for n in range(len(mz_vec_corr)-1):
      i = n
      if HS_num[n] == 0 and n != 0:
        HS_num[n] = n
        while Mod_HS_Dataframe['mod'][i+1]-Mod_HS_Dataframe['mod'][n] < hs_tol and exit_loop == 0:
          HS_num[i+1] = n
          HS_num[n] = n
          if i < (len(mz_vec_corr)-2):
            i = i+1
          else:
            exit_loop = 1

    if HS_num[-1] == 0:
      HS_num[-1] = n+1

    Mod_HS_Dataframe['HS Number'] = HS_num

    # Calculate number of members in HS
    HS_num_temp = np.unique(HS_num, return_counts = True)
    hsnumber = []

    for n in range(len(HS_num_temp[1])):
      for s in range(HS_num_temp[1][n]):
        hsnumber.append(HS_num_temp[1][n])

    Mod_HS_Dataframe['Homologues'] = hsnumber

    # Calculate corrected number of members in HS (same masses are neglected)
    HS_num2 = np.zeros(len(mz_vec_corr))
    for n in range(len(HS_num2)):
      HS_num2[n] = len(np.unique(round(Mod_HS_Dataframe['mz'][Mod_HS_Dataframe['HS Number'] == Mod_HS_Dataframe['HS Number'][n]])))

    Mod_HS_Dataframe['Unique Homologues'] = HS_num2

    # Check if number of homologues are greater than specified value
    HS_min_check = []
    for n in range(len(hsnumber)):
      if Mod_HS_Dataframe['Unique Homologues'][n] >= n_min:
        HS_min_check.append(True)
      else:
        HS_min_check.append(False)

    Mod_HS_Dataframe['min Homologues'] = HS_min_check

    # Check the features that were also detected by their fragments
    Mod_HS_Dataframe.sort_index(inplace=True)

    # Create KMD plots and output files
    Mod_HS_Dataframe.sort_values('min Homologues', inplace=True) #nach aufsteigendem modulo sortieren

    # Plot colored scatters & calculate normalized RT: Confirmed by HS
    Mod_HS_Dataframe_pos = Mod_HS_Dataframe[(Mod_HS_Dataframe['min Homologues'])]

    HS_unique = np.unique(Mod_HS_Dataframe_pos['HS Number'])
    RT_norm_fin = pd.Series(dtype='float')
    for i in range(len(HS_unique)):
      HS_temp = Mod_HS_Dataframe_pos[Mod_HS_Dataframe_pos['HS Number'] == HS_unique[i]]
      RT_min = np.amin(HS_temp['RT'])
      RT_max = np.amax(HS_temp['RT'])
      RT_norm = ((HS_temp['RT'] - RT_min) / (RT_max - RT_min))
      RT_norm_fin = pd.concat([RT_norm_fin,RT_norm])
    Mod_HS_Dataframe_pos = pd.concat([Mod_HS_Dataframe_pos,RT_norm_fin.rename('RT_norm3')],axis=1)

    # Write HS output file
    Mod_HS_Dataframe_pos_cp = Mod_HS_Dataframe_pos.copy()
    Mod_HS_Dataframe_pos_sort = Mod_HS_Dataframe_pos.sort_values(by='HS Number')
    Mod_HS_Dataframe_pos_sort.to_excel(os.path.join(Res_folder, frags_str[j] + '_homologous_series.xlsx'))

    # append all m/z's found by HS
    mz_HS_tot = np.append(mz_HS_tot,np.array(Mod_HS_Dataframe_pos_sort['mz']))

    # save only unique m/z's    
    mz_HS_tot = np.unique(mz_HS_tot)
    #=================================HSend======================================== 
    print('Job finished!\n')
    #globals().clear() # clear all variables to start the program new
    #==============================================================================    
