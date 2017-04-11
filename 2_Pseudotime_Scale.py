from __future__ import absolute_import

import random
from math import log
from operator import itemgetter

import numpy

import Orange
from Orange.bio import obiGEO
from Orange.bio.obiExpression import ExpressionSignificance_Test


# Utility functions

log2 = lambda x: log(x, 2.)

def my_ratio(x, y):
    """ compute the log-ratio """
    return log2(x/y)

def sign(x):
    return cmp(x, 0)

def common_domain(data1, data2):
    """Use only attributes that are in both data sets"""
    atts = sorted(set([a.name for a in data1.domain.attributes]).intersection(
               [a.name for a in data2.domain.attributes]))
    new1 = Orange.data.Table(Orange.data.Domain(atts + [data1.domain.classVar], data1.domain), data1)
    new2 = Orange.data.Table(Orange.data.Domain(atts + [data2.domain.classVar], data2.domain), data2)
    return new1, new2

def uniform_time_scale(attr_set):
    """ Obtains time points with a unique measure (the lowest among [min,h,d]) present in data"""
    ord_measures = ["min","h","d"]
    converter = {"min":{"h":60, "d":24*60}, "h": {"d":24}}
    measures = list(set([t.split(" ")[1] for t,s in attr_set]))
    if len(measures) == 1:
        time_points = [float(t.split(" ")[0]) for t,s in attr_set]
    else:
        first_measure = min([(ord_measures.index(m),m) for m in measures])[1]
        time_points = [float(t.split(" ")[0]) for t, s in attr_set
                       if t.split(" ")[1] == first_measure]
        time_points.extend([float(t.split(" ")[0]) * converter[first_measure][t.split(" ")[1]]
                            for t, s in attr_set if t.split(" ")[1] != first_measure])
    return time_points


def signed_PCA(data):
    pca = Orange.projection.linear.PCA(data, standardize=False, max_components=1)
    classifier = lambda X: [x[0].value for x in pca(X)]
    predictions = classifier(data)
    classes = [ex.getclass().value for ex in data]
    n = 0
    for i1,c1 in enumerate(classes):
        for i2,c2 in enumerate(classes[:i1]):
            n += cmp(c1,c2) * cmp(predictions[i1], predictions[i2])
    if n < 0:
        def invert(X):
            y = classifier(X)
            return -y if type(y) == float else [-x for x in y]
        spca = invert
    else:
        spca = classifier
    spca.pca = pca
    return spca

signed_PCA.name = 'PCA'


def conttime(data, d):
    """ Converts to continuous time labels """
    for a in data.domain.attributes:
        a.attributes['time'] = d[a.attributes['time']]

def conv(attr_set, ticks=True):
    """Obtain time points with a unique measure (the lowest among [min,h,d]) present in data"""
    ord_measures = ["min","h","d"]
    converter = {"min":{"h":60, "d":24*60}, "h": {"d":24}}
    measures = list(set([t.split(" ")[1] for t in attr_set]))
    if len(measures) == 1:
        time_points = [(t, float(t.split(" ")[0])) for t in attr_set]
    else:
        first_measure = min([(ord_measures.index(m),m) for m in measures])[1]
        time_points = [(t, float(t.split(" ")[0])) for t in attr_set if t.split(" ")[1] == first_measure]
        time_points.extend([(t, float(t.split(" ")[0]) * converter[first_measure][t.split(" ")[1]])
                            for t in attr_set if t.split(" ")[1] != first_measure])
    time_points.sort(key=itemgetter(1))
    if ticks:
        time_points = [(t[0],float(i)) for i,t in enumerate(time_points)]
    return dict(time_points)


def get_projections(data1, data2=None):
    labels1 = list(a.attributes['time'] for a in data1.domain.attributes)
    tdata1 = obiGEO.transpose(data1)
    if data2:
        labels2 = list('[%s]' % a.attributes['time'] for a in data2.domain.attributes)
        tdata2 = obiGEO.transpose(data2)
        tdata1, tdata2 = common_domain(tdata1, tdata2)
        classifier = signed_PCA(tdata1)
        proj1 = classifier(tdata1)
        proj2 = classifier(tdata2)
    else:
        classifier = signed_PCA(tdata1)
        proj1 = classifier(tdata1)
        proj2, labels2 = [], []
    return proj1, labels1, proj2, labels2, classifier


def scale_plot(names, points, filename='scale.png'):
    """Draw 'points' with 'names' for labels and save to 'filename'"""
    unames = set(names)
    ind = [sorted([i for i,n2 in enumerate(names) if n==n2], key=lambda x:points[x]) for n in sorted(unames)]
    ind.sort(key=lambda x:min(names[i] for i in x))

    medians = [numpy.median([points[ii[k]] for k in range(len(ii))]) for ii in ind]
    eps = float(max(points) - min(points)) / 75.0
    plt.clf()
    fig = plt.figure(figsize=(15,3), dpi=300)
    ax = fig.add_subplot(111)
    ax.arrow(min(points), -1, max(points)-min(points), 0, linewidth=1.2, head_width=0.5, head_length=1, fc='k', ec='k')
    height = 0
    rest = []
    ## (add colors if needed)
    colors = ['red', 'darkgreen', 'blue', 'darkorange', 'black']
    while ind!=[] or rest!=[]:
        if ind == []:
            ind, rest = rest, ind
            height += 1
        ax.plot([points[i] for i in ind[0]], [height]*len(ind[0]), color=colors[height], linewidth=1, marker='o', markersize=5)
        ax.plot([points[i] for i in ind[0]], [-3 for kk in range(len(ind[0]))], color=colors[height], linewidth=0, marker='|', markersize=5, mew=2)
        middle = ind[0][(len(ind[0])-1)/2]
        ##ax.text(points[middle], height+0.15, names[middle])
        lim = max(points[i] for i in ind[0])
        rest += [x for x in ind[1:] if min(points[i] for i in x) <= lim+eps]
        ind = [x for x in ind if min(points[i] for i in x) > lim+eps]
        ax.plot(medians[height], -1, color=colors[height], linewidth=0, marker='|', markersize=5, mew=1.2)
        ##ax.plot(medians[height], -1, color=colors[height], linewidth=0, marker='*', markersize=5, mew=1.2)
    plt.ylim(-3.5,height+1)
    plt.xlim(min(points)-3, max(points)+3)
    plt.yticks([])
    plt.setp(plt.gca(), 'xticklabels', [])
    ##plt.xlabel('Pseudotime')
    ##plt.grid(True)
    plt.savefig(filename, bbox_inches='tight')

def create_sampled_data_tr_ts(datax, at_list):
    """Sample model and test data using a bootstrap approach'"""
    namelabels = [at for at in datax.domain.attributes]
    labels1 = list(a.attributes['time'] for a in datax.domain.attributes)

    # for each stage, randomly sample the same number of samples with replacement
    # keep the remaining as a test set and save both data
    rand_all = []
    left_test = []
    for aa in at_list:
        namelab = [nl for i,nl in enumerate(namelabels) if labels1[i]==aa]
        numl = len(namelab)

        randlab = list(numpy.random.choice(namelab, numl, replace=True))
        ul = list(set(randlab))
        testsel = [atl for atl in namelab if atl not in randlab]
        rand_all.extend(randlab)
        left_test.extend(testsel)

    # save train data used to build the scale - leave sample names as they are
    data = []
    new_data = Orange.data.Table(datax.domain, data)
    for d in datax:
        ex = [d[at].value for at in rand_all]
        new_data.append(ex)
        new_data[-1]['gene'] = d['gene'].value

    # save test data - add data to a new Orange dataset
    # and use actual names of test samples to identify the collection time
    # and use it to measure the accuracy of sample placement on the scale
    lfnames = [atl.name for atl in left_test]
    h = [dd for dd in datax.domain if dd.name in lfnames]
    new_domaint = Orange.data.Domain(h, False)
    gg = Orange.feature.String('gene')
    new_domaint.add_meta(Orange.feature.Descriptor.new_meta_id(), gg)

    # create new dataset
    new_datat = Orange.data.Table(new_domaint, [])
    for d in datax:
        ex = [d[at].value for at in left_test]
        new_datat.append(ex)
        new_datat[-1]['gene'] = d['gene'].value

    return new_data, new_datat

# ------------------------------------------------------------------------------------------------------------------- #

if __name__ == '__main__':
    # 1) load data set 1
    data1 = Orange.data.Table('res/data_ALL_stages_batchcorr_svaseq_upper_quart_variant.tab')
    labels1 = list(a.attributes['time'] for a in data1.domain.attributes)
    attr_set = list(set(a.attributes['time'] for a in data1.domain.attributes))
    convd = conv(attr_set)
    conttime(data1, convd)

    # 2) project ona one dimensional line
    train = obiGEO.transpose(data1)
    classifier = signed_PCA(train)
    proj1 = classifier(train)
    classes = [str(ex.getclass().value) for ex in train]
    import matplotlib.pyplot as plt

    # 3) plot image and save projections
    scale_plot(classes, proj1, 'res/pseudot_projections_15_3_median_vline_nolab.pdf')

    # write cells projections
    h = [s['sample'].value for s in train]
    proj_cells = [(proj1[i],h[i]) for i in range(len(h))]
    proj_cells.sort()

    file_proj = open('res/cells_projections__.txt', 'w')
    file_proj.write('sample\tprojection\n')
    for pp in proj_cells:
        file_proj.write('%s\t%s\n' %(pp[1], pp[0]))
    file_proj.close()

    # 4) or, create a scale and project external data on the scale
    data1 = Orange.data.Table('res/data_ALL_stages_batchcorr_svaseq_upper_quart_variant.tab')
    data2 = Orange.data.Table('data/other_datasets/rpkm_3_month_svaseq_upper_quart_ourselgenes.tab')

    labels1 = list(a.attributes['time'] for a in data1.domain.attributes)
    attr_set = list(set(a.attributes['time'] for a in data1.domain.attributes))
    convd = conv(attr_set)
    conttime(data1, convd)
    train = obiGEO.transpose(data1)
    test = obiGEO.transpose(data2)
    cdtrain, cdtest = common_domain(train, test)

    classifier = signed_PCA(cdtrain)
    proj1 = classifier(cdtrain)
    proj2 = classifier(cdtest)
    h2 = [ssa.name for ssa in data2.domain.attributes]
    proj_cells2 = [(proj2[i],h2[i]) for i in range(len(h2))]
    proj_cells2.sort()

    ## proj2 from external data can be saved and plot as shown in 3)

    # Confidence interval computation: uses function create_sampled_data_tr_ts
    # to create 1000 bootstrap samples of the data, used for PCA model inference.
    # Remaining samples are projected as test sets and projections are saved
    # in order to be evaluated with accuracy measures.

    data1 = Orange.data.Table('res/data_ALL_stages_batchcorr_svaseq_upper_quart_variant.tab')
    nrep = 1000
    for krep in range(nrep):
        at_list = ['1 d', '7 d', '14 d', '21 d', '28 d']
        datay, dataz = create_sampled_data_tr_ts(data1, at_list)
        train = obiGEO.transpose(datay)
        test = obiGEO.transpose(dataz)
        classifier = signed_PCA(train)
        proj1 = classifier(train)
        classes = [str(ex.getclass().value) for ex in train]

        h = [s['sample'].value for s in train]
        proj_cells = [(proj1[i],h[i]) for i in range(len(h))]
        proj_cells.sort()

        proj2 = classifier(test)
        classest = [str(ex.getclass().value) for ex in test]

        h = [s['sample'].value for s in test]
        proj_cellst = [(proj2[i],h[i]) for i in range(len(h))]
        proj_cellst.sort()

        if (krep % 100) == 0:
            print krep
            file_proj = open('res/validation/cells_projections_test_bootstrap_train_test_%s.txt' %krep, 'w')

        file_proj.write('\t'.join([pp[1] for pp in proj_cellst])+ '\n')
        #del train, test, datay, dataz, proj_cellst
    file_proj.close()
