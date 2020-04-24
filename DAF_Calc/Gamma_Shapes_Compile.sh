# 4th Mar 2019
# Compiling gamma shapes into single file
# (to copy and paste in)

DAFS=($(seq 0.05 0.01 0.95))
NDAF=${#DAFS[@]}

rm gamma_shapes_HOL_highN0 gamma_shapes_HOL_lowN0
touch gamma_shapes_HOL_highN0 gamma_shapes_HOL_lowN0
for ((i=0; i<NDAF; i++));
do
cat HighNO/GSF${DAFS[${i}]}/sds_input.gamma_shapes >> gamma_shapes_HOL_highN0
cat LowNO/GSF${DAFS[${i}]}/sds_input.gamma_shapes >> gamma_shapes_HOL_lowN0
done
