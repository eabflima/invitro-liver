#!/bin/bash
rm *# *~ inp_dat.in &> /dev/null
time_t=2h
mito_c=10
folder_n=outputData
folder_e=(100 50 25 10 5)
for l in $(seq 1 5)
do
    echo "v_l_a_min = 0.0" | tee -a inp_dat.in &> /dev/null
    echo "v_l_a_max = 10.0" | tee -a inp_dat.in &> /dev/null
    echo "v_ic_min = 0.00001" | tee -a inp_dat.in &> /dev/null
    echo "v_ic_max = 1.0" | tee -a inp_dat.in &> /dev/null
    echo "v_std_min = 0.00001" | tee -a inp_dat.in &> /dev/null
    echo "v_std_max = 2.0" | tee -a inp_dat.in &> /dev/null
    echo "cal_name = ./all_data_trial9_${time_t}_${mito_c}.dat" | tee -a inp_dat.in &> /dev/null
    echo "n_samples = 4" | tee -a inp_dat.in &> /dev/null
    echo "n_cells = 5" | tee -a inp_dat.in &> /dev/null
    echo "cell_type = ${l}" | tee -a inp_dat.in &> /dev/null
    make run
    let n=l-1
    mv ./${folder_n} ./${folder_n}${folder_e[${n}]}
    rm *# *~ inp_dat.in &> /dev/null
done
mkdir Mito_${mito_c}
mv ./${folder_n}* ./Mito_${mito_c}/
mkdir ${time_t}
mv Mito_* ./${time_t}/
