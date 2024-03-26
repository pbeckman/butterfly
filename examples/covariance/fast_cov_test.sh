#!/bin/zsh
set -o nullglob

if [ "$#" -ne 1 ]; then
    echo "usage: mesh.obj"
    return
fi

MESH=$1
MESHBASE=${MESH##*/}
MESHNAME=${MESHBASE%.*}

# KAPPA, NU pairs to test
# PARAMS=(1e-4,0.0 1e-1,4.0 3e-1,0.5)
PARAMS=(3e-1,0.5)

NUM_SAMPLES=100

REFTOL=1e-4
# TOLS=(1e-1 1e-2 1e-3 $REFTOL)
TOLS=($REFTOL)
LOGPS=(2 4 6 8 10)

JULIA=true

# set up julia
if $JULIA; then julia --project=. -e "using Pkg; Pkg.instantiate()"; fi

# clean this directory and make output directory
rm -f -- *.bin *.txt
:> tols.txt
:> ps.txt
mkdir -p output

for logP in ${LOGPS[@]}
do
    P="$((2**$logP))"
    echo -n "$P " >> ps.txt
done

for TOL in ${TOLS[@]}
do
    echo -n "$TOL " >> tols.txt
done

for params in $PARAMS
do 
    IFS="," read KAPPA NU <<< "$params"

    echo "TOL\tRANK\tPRE(s)\t\tSAMP(s)\t\tCOMP(MB)\tUNCOMP(MB)\tUNTRUNC(MB)" >> $(printf "performance_kappa%.1e_nu%.1e.txt" $KAPPA $NU)

    # compute covariance for various butterfly tolerances
    for TOL in ${TOLS[@]}
    do
        # compute relevant quantities with butterfly covariance
        cmd="./lbo_cov $MESH $KAPPA $NU $NUM_SAMPLES $TOL" > $(printf "log_lbo_tol%.1e_kappa%.1e_nu%.1e.txt" $TOL $KAPPA $NU) 
        echo "\n$cmd"; eval $cmd

        # plot sample
        if $JULIA; then julia --project=. plot_sample.jl $MESH $KAPPA $NU lbo $TOL; fi
    done

    echo "\nP\tSAMP(s)" >> $(printf "performance_kappa%.1e_nu%.1e.txt" $KAPPA $NU)
    
    # compute covariances for various Chebyshev orders
    for logP in ${LOGPS[@]}
    do
        P="$((2**$logP))"

        # compute relevant quantities with Chebyshev covariance
        cmd="./cheb_cov $MESH $KAPPA $NU $NUM_SAMPLES $P" > $(printf "log_cheb_p%i_kappa%.1e_nu%.1e.txt" $P $KAPPA $NU)
        echo "\n$cmd"; eval $cmd

        # plot sample
        if $JULIA; then julia --project=. plot_sample.jl $MESH $KAPPA $NU cheb $P; fi
    done

    # plot error
    if $JULIA; then julia --project=. plot_fast_cov_test.jl $MESH $KAPPA $NU; fi
done

# # move generated files to output folder
# for file in *.bin *.txt
# do
#     mv $file output/"${MESHNAME}_$file"
# done
