# lameg_bstrm
Laminar MEG preprocessing and analysis tools implemented into the brainstorm environment. Original work for laminar MEG was implemented by Bonaiuto et al. 2018, 2021 in SPM. 

## Creating multilayer cortical meshes in brainstorm
Use `tess_create_multilayer_mesh` from the **anatomy** module to create a cortical mesh containing an arbitrary amount of layer. 
*Example Use:*
```
n_layers = 2;
wm_high_res = 'Path to *tess_cortex_white_high.m* in brainstorm_db';
wm_high_res = 'Path to *tess_cortex_pial_high.m* in brainstorm_db';
tess_create_multilayermesh(n_layers, wm_high_res, pial_high_res)
```

## Simulate single trial data 
Use `simulate_trial_data` from the **simulate** module to simulate n-Trials with activity on specific vertices.
*Example Use:*
```
data_struct = in_bst_data(*Path to Sensor File*);
headmodel_fname =*Path to headmodel*;
study_id = *ID to let brainstorm know where to save the data*; 

sim_params = struct( ...
                    'woi', [0.1 0.5],... %When the source is active
                    'foi', 20, ... % Simulated frequency (in Hz)
                    'sfreq', 250, ... % Sampling rate
                    'nTrials', 540, ... % Nr. of trials to simulate
                    'sim_loc', 7521, ... % Where to simulate on the cortical mesh
                    'snr_dB', -40); % Sensor-level signal to noise ratio

[sFiles, options] = simulate_trial_data(data_struct, headmodel_fname, sim_params, study_id);

```

## Extract Powerchanges to a baseline from all cortical sources
Use `extract_multilayer_power_changes.m` from the **inverse** module to extract a powerchange from a baseline to a window of interest on all cortical sources. Power is estimated as instanteneous power, calculated using `abs(hilbert(x))^2`.
*Example Use:*
```
sFiles = {*Dictionary of all sensor level files to be used*};
extract_multilayer_power_changes(sFiles);
```

## Compare Powerchanges between cortical layers
Use `compare_power_change_subject.m` from the **post_processing** module to compare previously extracted power changes between cortical layer.
*Example Use:*
```
sFiles = {*Dictionary of all sensor level files to be used*};

[contrast_name, t_vals, p_vals] = compare_power_change_subject(sFiles);
```
