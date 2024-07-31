DATA=$1

mkdir -p results_hpc

rsync -auv sesia@discovery.usc.edu:/home1/sesia/Workspace/collective-outlier-detection/code/third_party/nout/experiments/results/* results_hpc/
