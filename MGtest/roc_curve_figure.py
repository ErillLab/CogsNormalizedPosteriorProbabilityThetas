# Slightly modified from Pat O'Neill's repository code on github.
# Pat updated the computation of auc.

from math import sqrt,log,exp,pi,sin,cos,gamma,acos,sqrt
from matplotlib import pyplot as plt
import sys
import numpy as np
import itertools
def pairs(xs):
    return zip(xs[:-1],xs[1:])
def roc_curve(positives,negatives,color=None,annotate=False, plabel=None):
    """Given a set of positive scores and a set of negative scores, plot a
    ROC curve.

    Implements Algorithm 2 of:

    Fawcett, T. (2004), 'ROC graphs: Notes and practical
    considerations for researchers', ReCALL 31 (HPL-2003-4) , 1--38 .

    """
    instances = sorted([(x,0) for x in negatives] + [(x,1) for x in positives],
                       key = lambda (x,cls):x,
                       reverse=True)
    i = 0
    if color is None:
        color= 'b'
    tprs = []
    fprs = []
    thetas = []
    Np = float(len(positives))
    Nn = float(len(negatives))
    tp = 0
    fp = 0
    theta = min([x for (x,cls) in instances]) - 1
    theta_prev = theta
    for theta,cls in instances:
        #print "theta, cls", theta, cls
        if theta != theta_prev:
            tprs.append(tp/Np)
            fprs.append(fp/Nn)
            thetas.append(theta)
            theta_prev = theta
        if cls == 1:
            tp += 1
        else:
            fp += 1
    #print "tprs is", tp/Np, tp, Np
    #print "fprs is", fp/Nn, fp, Nn
    tprs.append(tp/Np)
    fprs.append(fp/Nn)
    thetas.append(theta)
    plt.plot(fprs,tprs,color=color,label=plabel)
    if annotate:
        theta_labels = ["%e" % theta for theta in thetas]
        annotations = unique(zip(fprs,tprs,theta_labels))
        modulus = len(annotations)/10
        print "%s unique annotations" % len(annotations)
        for i,(fpr,tpr,theta) in enumerate(annotations):
            if i % modulus == 0:
                plt.annotate(theta,xy=(fpr,tpr),xytext=(-20,20), textcoords = 'offset points',
                             ha = 'right', va = 'bottom',
                             bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
                             arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.xlabel("FPR", fontsize=22)
    plt.ylabel("TPR", fontsize=22)
    plt.plot([0,1],[0,1],linestyle='--')
    #plt.scatter(fprs, tprs)

    # add up squares below x0 and triangles above
    auc = sum((y0*(x1-x0)) + (0.5 *(y1-y0)*(x1-x0)) for ((x0,y0),(x1,y1)) in pairs(zip(fprs,tprs)))

    return fprs,tprs,thetas, auc
    
pscores1 = []
nscores1 = []
pscores2 = []
nscores2 = []
fname = sys.argv[1]
spos1 = int(sys.argv[2])
spos2 = int(sys.argv[3])
with open(fname, "r") as infile:
  inlines = infile.readlines()
for i in range(1, len(inlines)):
  line = inlines[i].strip()
  if len(line) < 10:
    continue
  items = line.split(",")
  if items[1] == "negative":
    nscores1.append(float(items[spos1]))
    nscores2.append(float(items[spos2]))
  else:
    pscores1.append(float(items[spos1]))
    pscores2.append(float(items[spos2]))
(fpr1, tpr1, th1, aucv1) = roc_curve(pscores1, nscores1, plabel="Permutation Test")
(fpr2, tpr2, th2, aucv2) = roc_curve(pscores2, nscores2, color='r', plabel="Bayesian Method")

plt.legend(loc='lower right')

plt.savefig("paper_figure3.png", dpi=300)
print "AUC LL:", aucv1
print "AUC PP:", aucv2
