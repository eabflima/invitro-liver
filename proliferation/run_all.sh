#!/bin/bash
rm *# *~ inp_dat.in &> /dev/null
folder_n=outputData
folder_e=(100 50 25 10 5 1)
for l in $(seq 1 5)
do
    let n=l-1
    echo -e "\t\t=====||${folder_n}${folder_e[${n}]}||====="
    apop_fm=`more ../apoptosis/2h/Mito_10/${folder_n}${folder_e[${n}]}/apop_mu_sigma.txt | awk NR==1`
    apop_fv=`more ../apoptosis/2h/Mito_10/${folder_n}${folder_e[${n}]}/apop_mu_sigma.txt | awk NR==2`
    ic_fm=`more ../apoptosis/2h/Mito_10/${folder_n}${folder_e[${n}]}/ic_mu_sigma.txt | awk NR==1`
    ic_fv=`more ../apoptosis/2h/Mito_10/${folder_n}${folder_e[${n}]}/ic_mu_sigma.txt | awk NR==2`
    rm *# *~ inp_dat.in &> /dev/null
    echo "v_l_p_min = 0.0" | tee -a inp_dat.in &> /dev/null
    echo "v_l_p_max = 10.0" | tee -a inp_dat.in &> /dev/null
    echo "v_l_k_min = 0.01" | tee -a inp_dat.in &> /dev/null
    echo "v_l_k_max = 1.0" | tee -a inp_dat.in &> /dev/null
    echo "v_std_min = 0.00001" | tee -a inp_dat.in &> /dev/null
    echo "v_std_max = 1.0" | tee -a inp_dat.in &> /dev/null
    echo "v_l_a_mean = ${apop_fm}" | tee -a inp_dat.in &> /dev/null
    echo "v_l_a_sigm = ${apop_fv}" | tee -a inp_dat.in &> /dev/null
    echo "v_ic_mean = ${ic_fm}" | tee -a inp_dat.in &> /dev/null
    echo "v_ic_sigm = ${ic_fv}" | tee -a inp_dat.in &> /dev/null
    echo "cal_name = ./all_data_trial8.dat" | tee -a inp_dat.in &> /dev/null
    echo "n_samples = 4" | tee -a inp_dat.in &> /dev/null
    echo "n_cells = 6" | tee -a inp_dat.in &> /dev/null
    echo "cell_type = ${l}" | tee -a inp_dat.in &> /dev/null
    make run
    mv ./${folder_n} ./${folder_n}${folder_e[${n}]}
done
mkdir Results
mv ./${folder_n}* ./Results/
