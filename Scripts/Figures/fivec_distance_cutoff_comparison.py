#!/usr/bin/env python

import sys
import os

import numpy
from pyx import canvas, text, path, graph, color, trafo, unit, attr, deco, style, bitmap

import hifive


unit.set(defaultunit="cm")
text.set(mode="latex")
text.preamble(r"\usepackage{times}")
text.preamble(r"\usepackage{sansmath}")
text.preamble(r"\sansmath")
text.preamble(r"\renewcommand*\familydefault{\sfdefault}")
painter = graph.axis.painter.regular( labeldist=0.1, labelattrs=[text.size(-3)], titleattrs=[text.size(-3)] )

methods = ['Prob', 'Exp', 'Bin', 'Exp-KR']
method_colors = {
        'Prob':color.cmyk.Black,
        'Exp':color.cmyk.CadetBlue,
        'Bin':color.cmyk.MidnightBlue, 
        'Exp-KR':color.cmyk.Mahogany,
}
meth_names = {"Prob":"Probability", "Exp":"Express", "Bin":"Binning",
              "Exp-KR":"ExpressKR"}


def main():
    out_fname = sys.argv[1]
    basedir = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-2])
    hic_fname1 = "%s/Data/HiC/HiFive/mm9_ESC_NcoI_prob.hcp" % basedir
    hic_fname2 = "%s/Data/HiC/HiFive/mm9_ESC_HindIII_prob.hcp" % basedir
    fivec_fnames = {
        "Prob_Phillips":"%s/Data/FiveC/HiFive/Phillips_ESC_prob.fcp" % basedir,
        "Prob_Nora":"%s/Data/FiveC/HiFive/Nora_ESC_male_E14_prob.fcp" % basedir,
        "Bin_Phillips":"%s/Data/FiveC/HiFive/Phillips_ESC_bin.fcp" % basedir,
        "Bin_Nora":"%s/Data/FiveC/HiFive/Nora_ESC_male_E14_bin.fcp" % basedir,
        "Exp_Phillips":"%s/Data/FiveC/HiFive/Phillips_ESC_exp.fcp" % basedir,
        "Exp_Nora":"%s/Data/FiveC/HiFive/Nora_ESC_male_E14_exp.fcp" % basedir,
        "Exp-KR_Phillips":"%s/Data/FiveC/HiFive/Phillips_ESC_expKR.fcp" % basedir,
        "Exp-KR_Nora":"%s/Data/FiveC/HiFive/Nora_ESC_male_E14_expKR.fcp" % basedir,
    }
    hic1 = hifive.HiC(hic_fname1)
    hic2 = hifive.HiC(hic_fname2)
    hic_hm = { 'Phillips': {} }
    fc = hifive.FiveC(fivec_fnames["Prob_Phillips"])
    fragments = fc.frags['fragments'][...]
    regions = fc.frags['regions'][...]
    for i in range(fc.frags['regions'].shape[0]):
            binbounds = numpy.hstack((
                    fragments['start'][regions['start_frag'][i]:regions['stop_frag'][i]].reshape(-1, 1),
                    fragments['stop'][regions['start_frag'][i]:regions['stop_frag'][i]].reshape(-1, 1)))
            binbounds = binbounds[numpy.where(fc.filter[regions['start_frag'][i]:regions['stop_frag'][i]])[0], :]
            hic_hm['Phillips'][i] = dynamically_bin(hic1, hic2, regions['chromosome'][i], binbounds)
    fc = hifive.FiveC(fivec_fnames["Prob_Nora"])
    fragments = fc.frags['fragments'][...]
    regions = fc.frags['regions'][...]
    binbounds = numpy.hstack((
            fragments['start'][regions['start_frag'][0]:regions['stop_frag'][0]].reshape(-1, 1),
            fragments['stop'][regions['start_frag'][0]:regions['stop_frag'][0]].reshape(-1, 1)))
    binbounds = binbounds[numpy.where(fc.filter[regions['start_frag'][0]:regions['stop_frag'][0]])[0], :]
    hic_hm['Nora'] = dynamically_bin(hic1, hic2, regions['chromosome'][0], binbounds)
    dist_corr = find_correlations( hic_hm, fivec_fnames, out_fname, True )
    nodist_corr = find_correlations( hic_hm, fivec_fnames, out_fname, False )
    c = canvas.canvas()
    width = 16.8
    spacer = 0.4
    plot_width = (width - spacer * 2) / 2.5
    plot_height = plot_width
    key_width = width - (plot_width + spacer) * 2 
    phillips_img = plot_correlation_diffs(dist_corr, nodist_corr, 'Phillips', plot_width, plot_height)
    nora_img = plot_correlation_diffs(dist_corr, nodist_corr, 'Nora', plot_width, plot_height)
    key_img = plot_key(key_width, plot_height)
    c.insert(phillips_img)
    c.insert(nora_img, [trafo.translate(plot_width + spacer, 0)])
    c.insert(key_img, [trafo.translate((plot_width + spacer) * 2, 0)])
    c.text(0, plot_height, "a",
       [text.halign.left, text.valign.top, text.size(-1)])
    c.text(plot_width + spacer, plot_height, "b",
       [text.halign.left, text.valign.top, text.size(-1)])
    c.writePDFfile(out_fname)

def find_correlations( hic_hm, fivec_fnames, out_fname, dist ):
    data = {}
    for meth in ['Prob', 'Bin', 'Exp', 'Exp-KR']:
        if dist:
            fname = fivec_fnames['%s_Phillips' % meth]
        else:
            fname = fivec_fnames['%s_Phillips' % meth].replace('.fcp', 'nodist.fcp')
        fc = hifive.FiveC(fname)
        fragments = fc.frags['fragments'][...]
        regions = fc.frags['regions'][...]
        counts = numpy.zeros(0, dtype=numpy.float64)
        expected = numpy.zeros(0, dtype=numpy.float64)
        hic_counts = numpy.zeros(0, dtype=numpy.float64)
        hic_expected = numpy.zeros(0, dtype=numpy.float64)
        skipped = []
        for i in range(fc.frags['regions'].shape[0]):
            temp = fc.cis_heatmap(i, datatype='fragment', arraytype='compact', binsize=0, skipfiltered=True)
            if temp is None:
                skipped.append(i)
                continue
            counts = numpy.hstack((counts, temp[:, :, 0].ravel()))
            expected = numpy.hstack((expected, temp[:, :, 1].ravel()))
            if meth == 'Prob':
                temp1 = numpy.zeros((temp.shape[0], temp.shape[1]), dtype=numpy.float32)
                temp1[numpy.where(temp[:, :, 0] > 0.0)] = 1
                binbounds = numpy.hstack((
                    fragments['start'][regions['start_frag'][i]:regions['stop_frag'][i]].reshape(-1, 1),
                    fragments['stop'][regions['start_frag'][i]:regions['stop_frag'][i]].reshape(-1, 1)))
                valid = numpy.where(fc.filter[regions['start_frag'][i]:regions['stop_frag'][i]])[0]
                binbounds = binbounds[valid, :]
                strands = fragments['strand'][regions['start_frag'][i]:regions['stop_frag'][i]][valid]
                temp = hic_hm['Phillips'][i][numpy.where(strands == 0)[0], :, :][:, numpy.where(strands == 1)[0], :]
                hic_counts = numpy.hstack((hic_counts, temp[:, :, 0].ravel()))
                hic_expected = numpy.hstack((hic_expected, temp[:, :, 1].ravel()))
        if meth == 'Prob':
            where = numpy.where(hic_expected > 0.0)[0]
            hic_counts[where] /= hic_expected[where]
            data["HiC_Phillips"] = numpy.copy(hic_counts)
        where = numpy.where(expected > 0.0)[0]
        counts[where] /= expected[where] 
        data["%s_Phillips" % meth] = numpy.copy(counts)
        if dist:
            fname = fivec_fnames['%s_Nora' % meth]
        else:
            fname = fivec_fnames['%s_Nora' % meth].replace('.fcp', 'nodist.fcp')
        fc = hifive.FiveC(fname)
        temp = fc.cis_heatmap(0, datatype='fragment', arraytype='compact', binsize=0, skipfiltered=True)
        counts = temp[:, :, 0].ravel()
        expected = temp[:, :, 1].ravel()
        if meth == 'Prob':
            temp1 = numpy.zeros((temp.shape[0], temp.shape[1]), dtype=numpy.float32)
            temp1[numpy.where(temp[:, :, 0] > 0.0)] = 1
            fragments = fc.frags['fragments'][...]
            regions = fc.frags['regions'][...]
            binbounds = numpy.hstack((
                    fragments['start'][regions['start_frag'][0]:regions['stop_frag'][0]].reshape(-1, 1),
                    fragments['stop'][regions['start_frag'][0]:regions['stop_frag'][0]].reshape(-1, 1)))
            binbounds = binbounds[numpy.where(fc.filter[regions['start_frag'][0]:regions['stop_frag'][0]])[0], :]
            strands = fragments['strand'][regions['start_frag'][0]:regions['stop_frag'][0]]
            temp = hic_hm['Nora'][numpy.where(strands==0)[0], :, :][:, numpy.where(strands == 1)[0], :]
            hic_counts = temp[:, :, 0].ravel()
            hic_expected = temp[:, :, 1].ravel()
            where = numpy.where(hic_expected > 0.0)[0]
            hic_counts[where] /= hic_expected[where]
            data["HiC_Nora"] = numpy.copy(hic_counts)
        where = numpy.where(expected > 0.0)[0]
        counts[where] /= expected[where] 
        data["%s_Nora" % meth] = numpy.copy(counts)
    correlations = {}
    for meth in methods:
        for name in ["Phillips", "Nora"]:
            valid = numpy.where((data["%s_%s" % (meth, name)] > 0.0) * (data["HiC_%s" % name] > 0.0))
            correlations["%s_%s" % (meth, name)] = numpy.corrcoef(numpy.log(data["%s_%s" % (meth, name)][valid]),
                                                              numpy.log(data["HiC_%s" % name][valid]))[0, 1]
    if not dist:
        output = open(out_fname.replace('pdf', 'txt'), 'w')
        print >> output, "Method\tPhillips\tNora"
        for meth in methods:
            temp = [meth]
            for name in ["Phillips", "Nora"]:
                temp.append(str(correlations["%s_%s" % (meth, name)]))
            print >> output, '\t'.join(temp)
        output.close()
    return correlations

def plot_correlation_diffs(corr1, corr2, name, width, height):
    ho = 1.2
    vo = 0.4
    plot_width = width - ho
    plot_height = height - vo
    diffs = {}
    ymin = numpy.inf
    ymax = -numpy.inf
    for n in ['Phillips', 'Nora']:
        for meth in meth_names.keys():
            cname = "%s_%s" % (meth, n)
            diff = corr2[cname] - corr1[cname]
            ymin = min(ymin, diff)
            ymax = max(ymax, diff)
    for meth in meth_names.keys():
        cname = "%s_%s" % (meth, name)
        diffs[meth] = corr2[cname] - corr1[cname]
    yspan = ymax - ymin
    ymin -= yspan * 0.05
    ymax += yspan * 0.05
    yspan = ymax - ymin
    c = canvas.canvas()
    g = graph.graphxy(width=plot_width, height=plot_height,
                      x=graph.axis.bar(painter=graph.axis.painter.bar(nameattrs=None)),
                      y=graph.axis.lin(painter=painter, min=ymin, max=ymax),
                      x2=graph.axis.lin(parter=None, min=0, max=1),
                      y2=graph.axis.lin(parter=None, min=0, max=1))
    w = plot_width / float(len(meth_names) + 1)
    y0 = -ymin / yspan * plot_height
    for i, meth in enumerate(methods):
        col = method_colors[meth]
        g.stroke( path.rect((i + 0.5) * w, y0, w, diffs[meth] / yspan * plot_height), [deco.filled([col])])
    g.stroke( path.line(0, y0, plot_width, y0), [style.linestyle.dotted, style.linewidth.THin])
    c.insert(g, [trafo.translate(ho, vo)])
    c.text(0, plot_height * 0.5 + vo, r"$r_{0K} - r_{50K}$",
           [text.halign.center, text.valign.top, text.size(-3), trafo.rotate(90)])
    c.text(plot_width * 0.5 + ho, vo * 0.5, name,
           [text.halign.center, text.valign.middle, text.size(-3)])
    return c

def plot_key(width, height):
    c = canvas.canvas()
    w = height / 7.0
    for i, meth in enumerate(methods):
        c.fill(path.rect(1.0, (6 - i) * w - 0.1, 0.2, 0.2), [method_colors[meth]])
        c.text(1.3, (6 - i) * w, "%s" % meth_names[meth],
               [text.halign.left, text.valign.middle, text.size(-3)])
    return c

def dynamically_bin(hic1, hic2, chrom, binbounds):
    unbinned1, map1 = hic1.cis_heatmap(chrom, start=binbounds[0, 0], stop=binbounds[-1, 1], datatype='fend',
                                 arraytype='full', returnmapping=True)
    unbinned2, map2 = hic2.cis_heatmap(chrom, start=binbounds[0, 0], stop=binbounds[-1, 1], datatype='fend',
                                 arraytype='full', returnmapping=True)
    map1[:, 2] = (map1[:, 0] + map1[:, 1])
    map2[:, 2] = (map2[:, 0] + map2[:, 1])
    allmap = numpy.vstack((map1, map2))
    allmap = allmap[numpy.argsort(allmap[:, 2]), :]
    indices1 = numpy.searchsorted(allmap[:, 2], map1[:, 2])
    indices1_1 = (indices1.reshape(-1, 1) * allmap.shape[0] + indices1.reshape(1, -1)).ravel()
    indices2 = numpy.searchsorted(allmap[:, 2], map2[:, 2])
    indices2_1 = (indices2.reshape(-1, 1) * allmap.shape[0] + indices2.reshape(1, -1)).ravel()
    unbinned = numpy.zeros((allmap.shape[0], allmap.shape[0], 2), dtype=numpy.float32)
    unbinned[:, :, 0] += numpy.bincount(indices1_1, minlength=allmap.shape[0] ** 2,
                                        weights=unbinned1[:, :, 0].ravel()).reshape(allmap.shape[0], -1)
    unbinned[:, :, 1] += numpy.bincount(indices1_1, minlength=allmap.shape[0] ** 2,
                                        weights=unbinned1[:, :, 1].ravel()).reshape(allmap.shape[0], -1)
    unbinned[:, :, 0] += numpy.bincount(indices2_1, minlength=allmap.shape[0] ** 2,
                                        weights=unbinned2[:, :, 0].ravel()).reshape(allmap.shape[0], -1)
    unbinned[:, :, 1] += numpy.bincount(indices2_1, minlength=allmap.shape[0] ** 2,
                                        weights=unbinned2[:, :, 1].ravel()).reshape(allmap.shape[0], -1)
    indices = numpy.triu_indices(allmap.shape[0], 1)
    unbinned = unbinned[indices[0], indices[1], :]
    binned, binmap = hic1.cis_heatmap(chrom, binbounds=binbounds, datatype='fend', arraytype='full',
                                      returnmapping=True)
    binned += hic2.cis_heatmap(chrom, binbounds=binbounds, datatype='fend', arraytype='full')
    indices = numpy.triu_indices(binbounds.shape[0], 1)
    upper = binned[indices[0], indices[1], :]
    hifive.hic_binning.dynamically_bin_cis_array(unbinned, allmap, upper, binmap,
                                                 expansion_binsize=0, minobservations=25)
    binned[indices[0], indices[1], :] = upper
    binned[indices[1], indices[0], :] = upper
    return binned

if __name__ == "__main__":
    main()
