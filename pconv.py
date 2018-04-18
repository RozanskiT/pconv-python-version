#!/usr/bin/python3
# -*- coding: utf-8 -*-
import sys, argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

"""
Simplified Python version of IDL program PCONV.pro of Ivan Hubeny & Thierry Lanz
http://nova.astro.umd.edu/
https://arxiv.org/abs/1706.01937
"""
def readModel(filename):
    model=pd.DataFrame()
    deepPoints=[]
    NDEPTH=0;NUMPAR=0;
    with open(filename) as f:
        line = f.readline()
        [NDEPTH,NUMPAR]=[int(s) for s in line.split()]
        while len(deepPoints) < NDEPTH:
            line = f.readline()
            l=[float(s.replace('D', 'E')) for s in line.split()]
            deepPoints.extend(l)
        dfSpectrum = pd.read_table(f,header=None,delim_whitespace=True)
        dfSpectrum = dfSpectrum.applymap(lambda x: float(x.replace('D', 'E')))
        model=np.transpose(dfSpectrum.values.flatten().reshape((-1,NUMPAR)))
        #model contains
        # in model[0,:] - Temperature
        # in model[1,:] - Electron density
        # in model[1,:] - Mass density
        # in the rest populations (unless it is model of acretion disc)
    return (NDEPTH,NUMPAR,deepPoints,model)

def main():
    descStr="""
    Display convergence logs of TLUSTY run.
    Needs files:
    *.7
    *.69
    *.9
    """

    parser=argparse.ArgumentParser(description=descStr)
    parser.add_argument(nargs='?',dest='modelName',default="fort",help="Core of model name, ex. fort; Default = 'fort' ")

    args=parser.parse_args()

    #read files
    dfCovergenceLog=pd.read_table(args.modelName+".9",header=0, skiprows=[0,2],delim_whitespace=True)
    dfTimingLog=pd.read_table(args.modelName+".69",header=None,delim_whitespace=True,)
    NDEPTH,NUMPAR,deepPoints,_=readModel(args.modelName+".7")
    ############################################################################
    #Process data a bit
    niter=max(dfCovergenceLog["ITER"])
    dfCovergenceLog=dfCovergenceLog.sort_values(by=["ITER","ID"])
    dfCovergenceLog["TEMP"]=dfCovergenceLog["TEMP"].apply(lambda x: float(x.replace('D', 'E')))
    dfCovergenceLog["MAXIMUM"]=dfCovergenceLog["MAXIMUM"].apply(lambda x: float(x.replace('D', 'E')))
    ############################################################################
    f, axarr = plt.subplots(2, 3)
    #PLOT TEMPERATURE
    for i in range(1,niter-1):
        axarr[0, 0].semilogx(deepPoints, dfCovergenceLog[dfCovergenceLog["ITER"]==i]["TEMP"])
    axarr[0, 0].set_title('Temperature')
    axarr[0, 0].set_ylabel('Relative Change')
    axarr[0, 0].set_xlabel('Column mass [g cm-2]')
    #END PLOT TEMPERATURE

    #PLOT TEMPERATURE LOGLOG
    for i in range(1,niter-1):
        axarr[1, 0].loglog(deepPoints, abs(dfCovergenceLog[dfCovergenceLog["ITER"]==i]["TEMP"]))
    axarr[1, 0].set_title('Temperature')
    axarr[1, 0].set_ylabel('Relative Change')
    axarr[1, 0].set_xlabel('Column mass [g cm-2]')
    #END PLOT TEMPERATURE LOGLOG

    #PLOT Maximum in State Vector
    for i in range(1,niter-1):
        axarr[0, 1].semilogx(deepPoints, dfCovergenceLog[dfCovergenceLog["ITER"]==i]["MAXIMUM"] )
    axarr[0, 1].set_title('Maximum in State Vector')
    axarr[0, 1].set_ylabel('Relative Change')
    axarr[0, 1].set_xlabel('Column mass [g cm-2]')
    #END PLOT  Maximum in State Vector

    #PLOT Maximum in State Vector
    for i in range(1,niter-1):
        axarr[1, 1].loglog(deepPoints, abs(dfCovergenceLog[dfCovergenceLog["ITER"]==i]["MAXIMUM"] ))
    axarr[1, 1].set_title('Maximum in State Vector')
    axarr[1, 1].set_ylabel('Relative Change')
    axarr[1, 1].set_xlabel('Column mass [g cm-2]')
    #END PLOT  Maximum in State Vector

    #PLOT Maximum of log(TEMP) reative change in each step
    axarr[0, 2].semilogy(np.arange(1,niter+1),[max(abs(dfCovergenceLog[dfCovergenceLog["ITER"]==i]["TEMP"].values)) for i in np.arange(1,niter+1)])
    axarr[0, 2].set_title('Temperature')
    axarr[0, 2].set_ylabel('Relative Change')
    axarr[0, 2].set_xlabel('Iteration')
    #END PLOT Maximum of log(TEMP) reative change in each step

    #PLOT  Maximum of log(MAXIMUM) reative change in each step
    axarr[1, 2].semilogy(np.arange(1,niter+1),[max(abs(dfCovergenceLog[dfCovergenceLog["ITER"]==i]["MAXIMUM"].values)) for i in np.arange(1,niter+1)])
    axarr[1, 2].set_title('Maximum in State Vector')
    axarr[1, 2].set_ylabel('Relative Change')
    axarr[1, 2].set_xlabel('Iteration')
    #END PLOT  Maximum of log(MAXIMUM) reative change in each step
    plt.show()

if __name__ == '__main__':
    main()
