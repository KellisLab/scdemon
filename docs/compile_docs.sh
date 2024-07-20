# conda activate scdtest
cd ..
pip install .

cd docs
make clean && make html
