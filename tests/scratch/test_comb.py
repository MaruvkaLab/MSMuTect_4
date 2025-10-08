import time, math

from src.IndelCalling.FisherTest import Fisher


def fisher_choose(iters: int):
    fisher = Fisher()
    for _ in range(iters):
        v = fisher.choose(1000, 79)


def inbuilt_choose(iters: int):
    fisher = Fisher()
    for _ in range(iters):
        v = math.comb(1000, 79)


def main():
    iters = int(1e6)
    st = time.time()
    inbuilt_choose(iters)
    e = time.time()
    print(e-st)

if __name__ == '__main__':
    main()