rm ./src/opca.o ./src/opca_binary.o ./src/opca_bimodal.o ./src/opca_pca.o ./lib/libopca.a
cd ./src

icpc -c opca_bimodal.cpp opca_binary.cpp opca_pca.cpp opca.cpp -I ../include -L../lib -llapack -lblas -lf2c -ltmg -lnlopt -lm
ar -cr libopca.a opca_bimodal.o opca_binary.o opca_pca.o opca.o
mv libopca.a ../lib
cd ..
icpc -o main_bimodal main_bimodal.cpp -L./lib -lopca -llapack -lblas -lf2c -ltmg -lnlopt -lm
icpc -o main_pca main_pca.cpp -L./lib -lopca -llapack -lblas -lf2c -lm

