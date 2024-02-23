Descriptions of .mat file data structure/organization are listed below for Figs. 2-9 of Aguirre, Woo et al., 2024 for the purpose of uploading to G-node.

Note: Select “AO_resample” or “SO_resample” substructures within each mat file for either action- or stimulus-based data, respectively.
Figure 2 – BLA inhibition impairs action-based probabilistic reversal learning, whereas vlOFC inhibition only slows early adjustment to reversals

For instructions for Figure 2-4,7 refer to: https://github.com/izquierdolab/trial-by-trial-pipeline.git 

Data files:
•	Figure 2C (Action-Based Discrimination): 
o	OFC hM4Di & EGFP learning curves: 'OFC_DIS_AB_OUT_DATA.mat'
•	AO_resample  tbl_trials 
	Virus==1 (OFC hM4Di) or Virus==0 (eGFP) 
o	BLA hM4Di learning curves: 'BLA_DIS_AB_OUT_DATA.mat'
	AO_resample  tbl_trials 
	Virus==1 (BLA hM4Di) 
•	Figure 2D (Action-Based Reversals 1 & 3): 
o	OFC hM4Di & EGFP learning curves: 'OFC_ALL_OUT_DATA.mat'
	AO_resample  tbl_trials 
	Virus==1 (OFC hM4Di) or Virus==0 (eGFP) 
	Drug=0 (VEH), =1 (CNO)
	Phase = 1 (Reversal 1), Phase =3 (Reversal 3) 
o	BLA hM4Di learning curves: 'BLA_ALL_OUT_DATA.mat'
	Virus==1 (BLA hM4Di) 
	Drug=0 (VEH), =1 (CNO)
	Phase = 1 (Reversal 1), Phase =3 (Reversal 3) 
Figure 3 - Female adjustment to reversals in the early phase was more affected by vlOFC inhibition than in males, whereas BLA’s role in probabilistic reversal learning was not sex-dependent

Data files:
•	Figure 3A (Action-Based Reversal 1): 
o	OFC hM4Di & EGFP by sex learning curves: 'OFC_ALL_OUT_DATA.mat'
	AO_resample  tbl_trials 
	Virus==1 (OFC hM4Di) or Virus==0 (eGFP) 
	Drug=0 (VEH), =1 (CNO)
	Phase = 1 (Reversal 1)
	Sex==1 (Male), Sex==0 (Female) 
o	BLA hM4Di by sex learning curves: 'BLA_ALL_OUT_DATA.mat'
	AO_resample  tbl_trials 
	Virus==1 (BLA hM4Di) 
	Drug=0 (VEH), =1 (CNO)
	Phase = 1 (Reversal 1)
	Sex==1 (Male), Sex==0 (Female) 
•	Figure 3B (Action-Based Reversal 3): 
o	OFC hM4Di & EGFP by sex learning curves: 'OFC_ALL_OUT_DATA.mat'
	AO_resample  tbl_trials 
	Virus==1 (OFC hM4Di) or Virus==0 (eGFP) 
	Drug=0 (VEH), =1 (CNO)
	Phase = 3 (Reversal 3)
	Sex==1 (Male), Sex==0 (Female) 
o	BLA hM4Di by sex learning curves: 'BLA_ALL_OUT_DATA.mat'
	AO_resample  tbl_trials 
	Virus==1 (BLA hM4Di) 
	Drug=0 (VEH), =1 (CNO)
	Phase = 1 (Reversal 1)
	Sex==1 (Male), Sex==0 (Female) 

Figure 4 – BLA inhibition impairs action-based probabilistic reversal learning in stimulus-naive animals, but does not affect initial action-based probabilistic discrimination or retention

Data files:
•	Figure 4A (Action-Based Reversal 3): 
o	BLA hM4Di learning curves: ‘BLA_AB1st_OUT_DATA.mat'
	AO_resample  tbl_trials 
	Virus==1 (BLA hM4Di) 
	Drug=0 (VEH), =1 (CNO)
	Phase==3 (Reversal 3 – 90/10)
o	BLA hM4Di by sex learning curves: ‘BLA_AB1st_OUT_DATA.mat'
	AO_resample  tbl_trials 
	Virus==1 (BLA hM4Di) 
	Drug=0 (VEH), =1 (CNO)
	Phase==3 (Reversal 3 – 90/10)
	Sex==1 (Male), Sex==0 (Females)
•	*Figure 4B (Action-Based Discrimination & Retention 90/10): 
o	BLA hM4Di learning curves: ‘BLA_AB1st_OUT_DATA.mat'

Figure 5 – The RL model fit to choice behavior indicates differential effects on inverse temperature and decay parameters between male and female rats during deterministic and probabilistic action-based reversals

Data files:
•	Figure 5A & 5B (Action-Based Reversal 1 & 3): 
o	 ‘output/model/fitNew/AO_IncludeNoChoice/RL1_decay.mat
	Use the plotting function ‘fit_functions/plot_figure_3_1.m’, section named %%Alternative version1.
	Set include_NoChoiceTrials = 1, since the model fitting include trials where choice is omitted.
	To run model fitting, run the fitting function fit_functions/run_fit.m. This uses the preprocessed behavioral data ‘allAnimalsDat_IncludeNoChoiceTrials.mat.'
Figure 6 – vlOFC, but not BLA, inhibition impairs adjustment to stimulus-based reversals as measured by probability correct

Data files:
•	Figure 6C (Stimulus-Based Discrimination): 
o	‘output/plot_data/block_output.mat’
	Use the plotting function ‘draw_all_plots.m,’ under the section with the function ‘figure_first_five_days.m.’
	Under %%settings, should select  taskType = "SO" in order to plot stimulus-outcome task. 
	Set MET.name = "pbetter" in order to plot for P(Correct).
•	Figure 6D (Stimulus-Based Reversal Adjustments/Transitions): 
o	‘output/plot_data/block_output.mat’
	Use the plotting function ‘draw_all_plots.m,’ under the section with the function ‘figure_days_around_rev.m.’
	Under %%settings, should select  taskType = "SO" in order to plot stimulus-outcome task. 
	Set MET.name = "pbetter" in order to plot for P(Correct).
o	To generate the above data file, run ‘produce_output_stats.m.’ This function uses the preprocessed behavioral data ‘allAnimalsDat.mat’ which oblivious to the no-choice trials. 
Figure 7 – BLA inhibition further slows incremental stimulus-based reversal learning whereas vlOFC inhibition abolishes learning after first reversal. Accuracy in stimulus learners measured by mean probability correct for the first five sessions of each deterministic (100/0) and probabilistic (90/10) reversal

Data files:
•	Figure 7 (Stimulus-Based Reversals 1-4 ): 
o	OFC hM4Di & EGFP learning curves: 'OFC_ALL_OUT_DATA.mat'
	SO_resample  tbl_day 
	Virus==1 (OFC hM4Di) or Virus==0 (eGFP) 
	Drug=0 (VEH), =1 (CNO)
	Phase = 1 (Reversal 1), =2 (Reversal 2), =3 (Reversal 3), =4 (Reversal 4)
o	BLA hM4Di learning curves: 'BLA_ALL_OUT_DATA.mat'
	SO_resample  tbl_day 
	Virus==1 (BLA hM4Di) 
	Drug=0 (VEH), =1 (CNO)
	Phase = 1 (Reversal 1), =2 (Reversal 2), =3 (Reversal 3), =4 (Reversal 4)
Figure 8 – BLA, but not vlOFC, inhibition significantly slowed action-based probabilistic reversal adjustment reflected in accuracy around reversals

Data files:
•	Figure 8 (Action-Based Reversal Adjustments/Transitions): 
o	‘output/plot_data/trial_output100.mat’
	Use the plotting function ‘draw_all_plots.m,’ under the section with the function ‘figure_trials_around_rev.m.’
	Under %%settings, should select  taskType = "AO" in order to plot action-outcome task. 
	Set MET.name = "pbetter" in order to plot for P(Correct).
o	To generate the data file, run ‘produce_output_stats.m.’ This function uses the preprocessed behavioral data ‘allAnimalsDat.mat’ which oblivious to the no-choice trials. 

Figure 9 – vlOFC, but not BLA, inhibition impairs adjustments to deterministic stimulus-based reversals as measured by changes in win–stay strategies

Data files:
•	Figure 9A (Action-Based Reversal Adjustments/Transitions: Win-Stay): 
o	‘output/plot_data/block_output.mat’
	Use the plotting function ‘draw_all_plots.m,’ under the section with the function ‘figure_days_around_rev.m.’
	Under %%settings, should select  taskType = "SO" in order to plot stimulus-outcome task. 
	Set MET.name = " winstay" in order to plot for Win-Stay.
•	Figure 9B (Action-Based Reversal Adjustments/Transitions: Lose-Shift): 
o	‘output/plot_data/block_output.mat’
	Set MET.name = " loseswitch" in order to plot for Lose-Shift, with the same instruction as above (Win-Stay).
o	To generate the data file, run ‘produce_output_stats.m.’ This function uses the preprocessed behavioral data ‘allAnimalsDat.mat’ which oblivious to the no-choice trials. 
