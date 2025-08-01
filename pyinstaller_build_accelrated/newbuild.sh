


# first build cython with bash build.sh
export PYTHONPATH=/home/avraham/MaruvkaLab/MSMuTect_0.5/
pyinstaller --onefile --name msmutect --hidden-import=tempfile --hidden-import=json --hidden-import=ctypes --hidden-import=platform \
    -p=/home/avraham/MaruvkaLab/MSMuTect_0.5/ \
    --collect-submodules=src.GenomicUtils --collect-submodules=src.Entry --collect-submodules=src.IndelCalling \
   /home/avraham/MaruvkaLab/MSMuTect_0.5/src/Entry/main.py
