
with open("/home/avraham/MaruvkaLab/Texas/strict_msmutect/results/res_low_purity.full.mut.tsv", 'r') as croc:
    fd = 1000
    found = 0
    while True:
        v=croc.readline()
        split_line = v.split("\t")
        call = split_line[50]
        if call.strip() == "GV":
            found+=1
            if len(split_line[3].strip())!=1:
                print(split_line[18])
                print(split_line[40])
                print(v)
                print("*************************************************")
            if found==fd:
                exit()


