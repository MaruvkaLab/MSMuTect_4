import pstats
from pstats import SortKey
# p = pstats.Stats('/home/avraham/MaruvkaLab/Texas/profiling_msmutect/post_removing_ifs2.prof')
# p = pstats.Stats('/home/avraham/MaruvkaLab/Texas/strict_msmutect/orig_msmutect.prof')
p = pstats.Stats("/home/avraham/MaruvkaLab/msmutect_runs/full_version/new_output.prof")
# p = pstats.Stats('/home/avraham/MaruvkaLab/Texas/strict_msmutect/strict_msmutect.prof')

p.strip_dirs().sort_stats(SortKey.CUMULATIVE).print_stats(50)
