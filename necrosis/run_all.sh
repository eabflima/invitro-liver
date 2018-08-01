#!/bin/bash
rm *# *~ inp_dat.in &> /dev/null
folder_n=outputData
folder_e=(100 50 25)
nut_c=(75 50 25 0)
pnut_c=(0.75 0.5 0.25 0.0)
#pnut_c=(1.0 1.0 1.0 1.0)
for c in $(seq 0 3)
do
    for l in $(seq 1 3)
    do
        let n=l-1
        echo -e "\t\t=====||${folder_n}${folder_e[${n}]}||====="
        apop_fm=`more ../proliferation/Results/${folder_n}${folder_e[${n}]}/apop_mu_sigma.txt | awk NR==1`
        apop_fv=`more ../proliferation/Results/${folder_n}${folder_e[${n}]}/apop_mu_sigma.txt | awk NR==2`
        ic_fm=`more ../proliferation/Results/${folder_n}${folder_e[${n}]}/ic_mu_sigma.txt | awk NR==1`
        ic_fv=`more ../proliferation/Results/${folder_n}${folder_e[${n}]}/ic_mu_sigma.txt | awk NR==2`
        prol_fm=`more ../proliferation/Results/${folder_n}${folder_e[${n}]}/prol_mu_sigma.txt | awk NR==1`
        prol_fv=`more ../proliferation/Results/${folder_n}${folder_e[${n}]}/prol_mu_sigma.txt | awk NR==2`
        supp_fm=`more ../proliferation/Results/${folder_n}${folder_e[${n}]}/supp_mu_sigma.txt | awk NR==1`
        supp_fv=`more ../proliferation/Results/${folder_n}${folder_e[${n}]}/supp_mu_sigma.txt | awk NR==2`
        rm *# *~ inp_dat.in &> /dev/null
        echo "v_l_n_min = 0.0" | tee -a inp_dat.in &> /dev/null
        echo "v_l_n_max = 10.0" | tee -a inp_dat.in &> /dev/null
        echo "v_std_min = 0.00001" | tee -a inp_dat.in &> /dev/null
        echo "v_std_max = 1.0" | tee -a inp_dat.in &> /dev/null
        echo "v_l_p_mean = ${prol_fm}" | tee -a inp_dat.in &> /dev/null
        echo "v_l_p_sigm = ${prol_fv}" | tee -a inp_dat.in &> /dev/null
        echo "v_l_k_mean = ${supp_fm}" | tee -a inp_dat.in &> /dev/null
        echo "v_l_k_sigm = ${supp_fv}" | tee -a inp_dat.in &> /dev/null
        echo "v_l_a_mean = ${apop_fm}" | tee -a inp_dat.in &> /dev/null
        echo "v_l_a_sigm = ${apop_fv}" | tee -a inp_dat.in &> /dev/null
        echo "v_ic_mean = ${ic_fm}" | tee -a inp_dat.in &> /dev/null
        echo "v_ic_sigm = ${ic_fv}" | tee -a inp_dat.in &> /dev/null
        echo "v_nut_val = ${pnut_c[${c}]}" | tee -a inp_dat.in &> /dev/null
        echo "cal_name = ./all_data_trial12_${nut_c[${c}]}.dat" | tee -a inp_dat.in &> /dev/null
        echo "n_samples = 4" | tee -a inp_dat.in &> /dev/null
        echo "n_cells = 3" | tee -a inp_dat.in &> /dev/null
        echo "cell_type = ${l}" | tee -a inp_dat.in &> /dev/null
        make run
        mv ./${folder_n} ./${folder_n}${folder_e[${n}]}
    done
    mkdir FBS_${nut_c[${c}]}
    mv ./${folder_n}* ./FBS_${nut_c[${c}]}/
done    
mkdir Results
mv ./FBS_* ./Results/
