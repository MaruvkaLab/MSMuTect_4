

#ReadsFetcher.cpython-38-x86_64-linux-gnu.so
#LocusFile.cpython-38-x86_64-linux-gnu.so
#CallAlleles.cpython-38-x86_64-linux-gnu.so
#Locus.cpython-38-x86_64-linux-gnu.so
#Histogram.cpython-38-x86_64-linux-gnu.so
#CallMutations.cpython-38-x86_64-linux-gnu.so
#AlleleSet.cpython-38-x86_64-linux-gnu.so
#FisherTest.cpython-38-x86_64-linux-gnu.so
#PairFileBatches.cpython-38-x86_64-linux-gnu.so
#SingleFileBatches.cpython-38-x86_64-linux-gnu.so

#CURRENT_PATH=$(readlink -f "$0")
#BUILD_DIR=$(dirname $CURRENT_PATH)
#MSMUTECT_DIR=$(dirname $BUILD_DIR)
#echo $MSMUTECT_DIR
export PYTHONPATH=/home/avraham/MaruvkaLab/MSMuTect_0.5/
pyinstaller --onefile --name msmutect --hidden-import=tempfile --hidden-import=json --hidden-import=ctypes --hidden-import=platform \
    --add-binary="/home/avraham/MaruvkaLab/MSMuTect_0.5/src/GenomicUtils/ReadsFetcher.cpython-38-x86_64-linux-gnu.so:src/GenomicUtils/msmutect/ReadsFetcher.cpython-38-x86_64-linux-gnu.so" \
    --add-binary="/home/avraham/MaruvkaLab/MSMuTect_0.5/src/GenomicUtils/LocusFile.cpython-38-x86_64-linux-gnu.so:src/GenomicUtils/LocusFile.cpython-38-x86_64-linux-gnu.so" \
    --add-binary="/home/avraham/MaruvkaLab/MSMuTect_0.5/src/IndelCalling/CallAlleles.cpython-38-x86_64-linux-gnu.so:src/IndelCalling/CallAlleles.cpython-38-x86_64-linux-gnu.so" \
    --add-binary="/home/avraham/MaruvkaLab/MSMuTect_0.5/src/IndelCalling/Locus.cpython-38-x86_64-linux-gnu.so:src/IndelCalling/Locus.cpython-38-x86_64-linux-gnu.so" \
    --add-binary="/home/avraham/MaruvkaLab/MSMuTect_0.5/src/IndelCalling/Histogram.cpython-38-x86_64-linux-gnu.so:src/IndelCalling/Histogram.cpython-38-x86_64-linux-gnu.so" \
    --add-binary="/home/avraham/MaruvkaLab/MSMuTect_0.5/src/IndelCalling/CallMutations.cpython-38-x86_64-linux-gnu.so:src/IndelCalling/CallMutations.cpython-38-x86_64-linux-gnu.so" \
    --add-binary="/home/avraham/MaruvkaLab/MSMuTect_0.5/src/IndelCalling/AlleleSet.cpython-38-x86_64-linux-gnu.so:src/IndelCalling/AlleleSet.cpython-38-x86_64-linux-gnu.so" \
    --add-binary="/home/avraham/MaruvkaLab/MSMuTect_0.5/src/IndelCalling/FisherTest.cpython-38-x86_64-linux-gnu.so:src/IndelCalling/FisherTest.cpython-38-x86_64-linux-gnu.so" \
    --add-binary="/home/avraham/MaruvkaLab/MSMuTect_0.5/src/Entry/PairFileBatches.cpython-38-x86_64-linux-gnu.so:src/Entry/PairFileBatches.cpython-38-x86_64-linux-gnu.so" \
    --add-binary="/home/avraham/MaruvkaLab/MSMuTect_0.5/src/Entry/SingleFileBatches.cpython-38-x86_64-linux-gnu.so:src/Entry/SingleFileBatches.cpython-38-x86_64-linux-gnu.so" \
    /home/avraham/MaruvkaLab/MSMuTect_0.5/src/Entry/main.py
