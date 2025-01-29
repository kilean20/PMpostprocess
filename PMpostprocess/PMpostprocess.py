from .signal_processing import process_profile_signal, project_L3_to_xyuv
from .profilemonitor_plot_util import profilemonitor_plot as pmplt
pmplt = pmplt()
import os
import re
import numpy as np
import pandas as pd


def list_raw_files(data_path):
    path = os.path.join(data_path, 'data')
    pattern = re.compile(r".*_PM_D\d{4}_\d{8}_\d{6}\.dat$")
    files = [f for f in os.listdir(path) if pattern.match(f)]
    return files


def read_raw_file(fname):
    # Retrieve parameters and raw data
    par_dict = pmplt.get_param(fname)
    raw = pmplt.get_rawdata(fname, par_dict)
    
    data = []
    
    for i in range(1, 4):  # Process for indices 1, 2, 3
        # Extract and sort positions and signals
        pos = np.array(raw[f'a1p{i}'].split()[4:], dtype=float)
        sig = np.array(raw[f'a1s{i}'].split()[4:], dtype=float)
        idx = np.argsort(pos)
        pos = pos[idx]
        sig = sig[idx]
        
        # Identify largest gaps in positions
        gaps = pos[1:] - pos[:-1]
        largest_gaps = np.sort(np.argsort(gaps)[-2:])
        
        # Find the domain containing the peak signal
        peak = np.argmax(sig)
        
        if largest_gaps[1] < peak:  # Peak is after the second largest gap
            pos = pos[largest_gaps[1] + 1:]
            sig = sig[largest_gaps[1] + 1:]
        elif largest_gaps[0] < peak:  # Peak is between the two largest gaps
            pos = pos[largest_gaps[0] + 1:largest_gaps[1]]
            sig = sig[largest_gaps[0] + 1:largest_gaps[1]]
        else:  # Peak is before the first largest gap
            pos = pos[:largest_gaps[0]]
            sig = sig[:largest_gaps[0]]
        # Store the filtered data
        data.append(np.stack((pos, sig)))
        
    par_dict['lraw'] = data
    return par_dict


def read_raw_data_from_directory(data_path):
    files = list_raw_files(data_path)
    raw = {}
    for f in files:
        par_dict = read_raw_file(os.path.join(data_path, 'data', f))
        raw[f] = par_dict   
    return raw
    
    
def read_all_raw_data(ISAAC_helper,ISAAC_data_rel_path,ISAAC_database_path,within_minutes=240):
    raw= read_raw_data_from_directory(os.path.join(ISAAC_database_path,ISAAC_data_rel_path))
    for rel_path in ISAAC_helper.get_related_ISAAC_data_rel_path(ISAAC_data_rel_path,ISAAC_database_path,within_minutes=within_minutes):
        raw.update(read_raw_data_from_directory(os.path.join(ISAAC_database_path,rel_path)))
    return raw
    

def reprocess_all_raw_data(raw):
    for key, value in raw.items():
        r = 0.5 * value['wire_diam']
        value['postprocess'] = {'denoised': [],
                                'projected': {}}
        rms = []
        arr_rms = []
        for ic, coord in enumerate(value['coord'][1:]):
            x = value['lraw'][ic][0]
            y = value['lraw'][ic][1]
            tmp = process_profile_signal(x, y, r)
            value['postprocess']['denoised'].append(tmp)
            rms.append(tmp['rms_deconv'])
            arr_rms.append(tmp['MC_stat']['arr_deconv_rms'])
        value.pop('lraw')
        
        x,y,cxy,u,v = project_L3_to_xyuv(*rms,value)
        arr_x,arr_y,arr_cxy,arr_u,arr_v = project_L3_to_xyuv(*arr_rms,value)
        value['postprocess']['projected']['xrms']=x
        value['postprocess']['projected']['yrms']=y
        value['postprocess']['projected']['cxy']=cxy
        value['postprocess']['projected']['urms']=u
        value['postprocess']['projected']['vrms']=v
        value['postprocess']['projected']['xrms_err']=np.std(arr_x)
        value['postprocess']['projected']['yrms_err']=np.std(arr_y)
        value['postprocess']['projected']['cxy_err']=np.std(arr_cxy)
        value['postprocess']['projected']['u_err']=np.std(arr_u)
        value['postprocess']['projected']['v_err']=np.std(arr_v)
    return raw


            
def reprocess_CS(FLAME_helper,ISAAC_data,raw_data):

    reconst_summary = ISAAC_data['reconst_summary']
    scan_data = reconst_summary['scan_data']
    nscan = len(scan_data)
    raw = raw_data
    raw_data_flist = [scan_data[i]['res_monitorl'][0]['file'][3:] for i in range(nscan)]
    raw = {key: raw[key] for key in raw.keys() & raw_data_flist}
    assert raw.keys() & raw_data_flist
    raw = reprocess_all_raw_data(raw)
    
    reconst_input = ISAAC_data['reconst_input']
    select = reconst_input['reconst_input']['measurement_select'][-1]
    
    fm_path = ISAAC_data['fmlatfile']
    fm_orig = FLAME_helper.ModelFlame(fm_path)
    fm = FLAME_helper.ModelFlame(fm_path)
    fm_evals = {}
    fm_goals = {}
    fm_goals_postprocessed = {}
    fm_goals_err = {}

    for i in range(nscan):
        for fm_elem,val_field in reconst_input['meas'][i]['flamevall'].items():
            k = (fm_elem,val_field[1])
            if k not in fm_evals:
                fm_evals[k] = [None]*nscan
            fm_evals[k][i] = val_field[0]
            
        for fm_elem,meas in reconst_input['meas'][i]['monitorl'].items():
            for goal in ['xrms','yrms','cxy']:
                k = (fm_elem,goal)
                if k not in fm_goals_postprocessed:
                    fm_goals[k] = [None]*nscan
                    fm_goals_postprocessed[k] = [None]*nscan
                    fm_goals_err[k] = [None]*nscan
                # if selected, record goal
                if select[goal][i]: 
                    fm_goals[k][i] = meas[goal]
                    fm_goals_postprocessed[k][i] = raw[raw_data_flist[i]]['postprocess']['projected'][goal]
                    fm_goals_err[k][i] = raw[raw_data_flist[i]]['postprocess']['projected'][goal+'_err']
                    
    fm_evals = FLAME_helper.make_FLAME_evals_or_goals(fm,df=pd.DataFrame(fm_evals))
    fm_goals = FLAME_helper.make_FLAME_evals_or_goals(fm,df=pd.DataFrame(fm_goals))
    fm_goals_postprocessed = FLAME_helper.make_FLAME_evals_or_goals(fm,df=pd.DataFrame(fm_goals_postprocessed))

    fm_goals_err = pd.DataFrame(fm_goals_err)
    goal_std = fm_goals_postprocessed['df'].std()
    goal_err_mean = fm_goals_err.mean()
    
    scaler = np.nanmean(goal_std/(goal_err_mean + 2e-4*goal_std))
    fm_goals_postprocessed['err'] = fm_goals_err
    fm_goals_postprocessed['normalization_factor'] = fm_goals_err*scaler  # make loss roughly order of 1
    
    FLAME_helper.fit_moment1(fm,fm_evals,fm_goals_postprocessed,stop_criteria=0.1)
    monitor_indices = [v['index'] for v in fm_goals_postprocessed['info'].values()]
    fm_eval_result = FLAME_helper.evaluate_flame_evals(fm_evals,fm,to_element=max(monitor_indices)+1, monitor_indices=monitor_indices)
    
    CS_new = FLAME_helper.bmstate2cs(fm.bmstate)
    CS_orig = FLAME_helper.bmstate2cs(fm_orig.bmstate)
    #CS_df = pd.DataFrame([CS_new, CS_orig], index=['CS_new', 'CS_orig'])
    
    return CS_new, CS_orig, fm, fm_evals, fm_goals, fm_goals_postprocessed, fm_eval_result, raw, raw_data_flist
    
    
