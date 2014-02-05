%PWMrandomizeBayes_run('PCM.txt',3,10);
%mcc PWMrandomizeBayes_run -W main -L Cpp -t -T link:exe -h -v -a ./ 
%mcc -T link:exe PWMrandomizeBayes_run_main.o PWMrandomizeBayes_run_mcc_component_data.o 
%mcc -m -v -R nojvm PWMrandomizeBayes_run -h -a ./
mcc -B sgl -v  PWMrandomizeBayes_run -a ./



