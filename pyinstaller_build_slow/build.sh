


# first build cython with bash build.sh
export PYTHONPATH=/home/avraham/MaruvkaLab/MSMuTect_0.5/
pyinstaller --onefile --name msmutect --hidden-import=tempfile --hidden-import=json --hidden-import=ctypes --hidden-import=platform \
    -p=/home/avraham/MaruvkaLab/MSMuTect_0.5/ \
    --collect-submodules=src.GenomicUtils --collect-submodules=src.Entry --collect-submodules=src.IndelCalling \
    --add-binary="/home/avraham/MaruvkaLab/MSMuTect_0.5/src/GenomicUtils/ReadsFetcher.cpython-38-x86_64-linux-gnu.so:src/GenomicUtils/" \
    --add-binary="/home/avraham/MaruvkaLab/MSMuTect_0.5/src/GenomicUtils/LocusFile.cpython-38-x86_64-linux-gnu.so:src/GenomicUtils/" \
    --add-binary="/home/avraham/MaruvkaLab/MSMuTect_0.5/src/IndelCalling/CallAlleles.cpython-38-x86_64-linux-gnu.so:src/IndelCalling/" \
    --add-binary="/home/avraham/MaruvkaLab/MSMuTect_0.5/src/IndelCalling/Locus.cpython-38-x86_64-linux-gnu.so:src/IndelCalling/" \
    --add-binary="/home/avraham/MaruvkaLab/MSMuTect_0.5/src/IndelCalling/Histogram.cpython-38-x86_64-linux-gnu.so:src/IndelCalling/" \
    --add-binary="/home/avraham/MaruvkaLab/MSMuTect_0.5/src/IndelCalling/CallMutations.cpython-38-x86_64-linux-gnu.so:src/IndelCalling/" \
    --add-binary="/home/avraham/MaruvkaLab/MSMuTect_0.5/src/IndelCalling/AlleleSet.cpython-38-x86_64-linux-gnu.so:src/IndelCalling/" \
    --add-binary="/home/avraham/MaruvkaLab/MSMuTect_0.5/src/IndelCalling/FisherTest.cpython-38-x86_64-linux-gnu.so:src/IndelCalling/" \
    --add-binary="/home/avraham/MaruvkaLab/MSMuTect_0.5/src/Entry/PairFileBatches.cpython-38-x86_64-linux-gnu.so:src/Entry/" \
    --add-binary="/home/avraham/MaruvkaLab/MSMuTect_0.5/src/Entry/SingleFileBatches.cpython-38-x86_64-linux-gnu.so:src/Entry/" \
    /home/avraham/MaruvkaLab/MSMuTect_0.5/src/Entry/main.py
