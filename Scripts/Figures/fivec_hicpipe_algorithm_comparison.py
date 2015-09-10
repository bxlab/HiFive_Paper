#!/usr/bin/env python

import sys
import os

import numpy
from pyx import canvas, text, path, graph, color, trafo, unit, attr, deco, style, bitmap
import h5py

import hifive


unit.set(defaultunit="cm")
text.set(mode="latex")
text.preamble(r"\usepackage{times}")
text.preamble(r"\usepackage{sansmath}")
text.preamble(r"\sansmath")
text.preamble(r"\renewcommand*\familydefault{\sfdefault}")
painter = graph.axis.painter.regular( labeldist=0.1, labelattrs=[text.size(-3)], titleattrs=[text.size(-3)] )

methods = ['Raw', 'Prob', 'Exp', 'Bin', 'Exp-KR']
method_colors = {
        'Prob':color.cmyk.Black,
        'Exp':color.cmyk.CadetBlue,
        'Bin':color.cmyk.MidnightBlue, 
        'Raw':color.cmyk.Dandelion,
        'Exp-KR':color.cmyk.Mahogany,
}


def main():
    out_fname = sys.argv[1]
    basedir = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-2])
    hic_phillips_fname1 = "%s/Data/HiC/HiCPipe/HM/mm9_ESC_NcoI_Phillips.hch" % basedir
    hic_phillips_fname2 = "%s/Data/HiC/HiCPipe/HM/mm9_ESC_HindIII_Phillips.hch" % basedir
    hic_nora_fname1 = "%s/Data/HiC/HiCPipe/HM/mm9_ESC_NcoI_Nora.hch" % basedir
    hic_nora_fname2 = "%s/Data/HiC/HiCPipe/HM/mm9_ESC_HindIII_Nora.hch" % basedir
    hic_phillips1 = h5py.File(hic_phillips_fname1, 'r')
    hic_phillips2 = h5py.File(hic_phillips_fname2, 'r')
    hic_nora1 = h5py.File(hic_nora_fname1, 'r')
    hic_nora2 = h5py.File(hic_nora_fname2, 'r')
    hm_phillips = {}
    hm_nora = {}
    for key in hic_phillips1.keys():
        if key.count('unbinned_counts') == 0:
            continue
        region = int(key.split('.')[0])
        hm_phillips[region] = dynamically_bin(hic_phillips1, hic_phillips2, region)
    for key in hic_nora1.keys():
        if key.count('unbinned_counts') == 0:
            continue
        region = int(key.split('.')[0])
        hm_nora[region] = dynamically_bin(hic_nora1, hic_nora2, region)

    fivec_fnames = {
        "Prob_Phillips":"%s/Data/FiveC/HiFive/Phillips_ESC_probnodist.fcp" % basedir,
        "Prob_Nora":"%s/Data/FiveC/HiFive/Nora_ESC_male_E14_probnodist.fcp" % basedir,
        "Bin_Phillips":"%s/Data/FiveC/HiFive/Phillips_ESC_binnodist.fcp" % basedir,
        "Bin_Nora":"%s/Data/FiveC/HiFive/Nora_ESC_male_E14_binnodist.fcp" % basedir,
        "Exp_Phillips":"%s/Data/FiveC/HiFive/Phillips_ESC_expnodist.fcp" % basedir,
        "Exp_Nora":"%s/Data/FiveC/HiFive/Nora_ESC_male_E14_expnodist.fcp" % basedir,
        "Exp-KR_Phillips":"%s/Data/FiveC/HiFive/Phillips_ESC_expKRnodist.fcp" % basedir,
        "Exp-KR_Nora":"%s/Data/FiveC/HiFive/Nora_ESC_male_E14_expKRnodist.fcp" % basedir,
    }
    data = {}
    imgs = {}
    ratio1 = 0
    ratio2 = 0
    for meth in ['Prob', 'Bin', 'Exp', 'Exp-KR']:
        fc = hifive.FiveC(fivec_fnames["%s_Phillips" % meth])
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
            if i == 6:
                ratio1 = temp.shape[1] / float(temp.shape[0])
                imgs["%s_Phillips" % meth] = hifive.plotting.plot_full_array(temp, symmetricscaling=False)
            if meth == 'Prob':
                temp1 = numpy.zeros((temp.shape[0], temp.shape[1]), dtype=numpy.float32)
                temp1[numpy.where(temp[:, :, 0] > 0.0)] = 1
                if i == 6:
                    imgs["Raw_Phillips"] = hifive.plotting.plot_full_array(
                            numpy.dstack((temp[:, :, 0], temp1)), symmetricscaling=False)
                binbounds = numpy.hstack((
                    fragments['start'][regions['start_frag'][i]:regions['stop_frag'][i]].reshape(-1, 1),
                    fragments['stop'][regions['start_frag'][i]:regions['stop_frag'][i]].reshape(-1, 1)))
                valid = numpy.where(fc.filter[regions['start_frag'][i]:regions['stop_frag'][i]])[0]
                binbounds = binbounds[valid, :]
                temp = hm_phillips[i]
                strands = fragments['strand'][regions['start_frag'][i]:regions['stop_frag'][i]][valid]
                temp = temp[numpy.where(strands == 0)[0], :, :][:, numpy.where(strands == 1)[0], :]
                hic_counts = numpy.hstack((hic_counts, temp[:, :, 0].ravel()))
                hic_expected = numpy.hstack((hic_expected, temp[:, :, 1].ravel()))
                if i == 6:
                    imgs["HiC_Phillips"] = hifive.plotting.plot_full_array(temp, symmetricscaling=False)
        if meth == 'Prob':
            data["Raw_Phillips"] = numpy.copy(counts)
            where = numpy.where(hic_expected > 0.0)[0]
            hic_counts[where] /= hic_expected[where]
            data["HiC_Phillips"] = numpy.copy(hic_counts)
        where = numpy.where(expected > 0.0)[0]
        counts[where] /= expected[where] 
        data["%s_Phillips" % meth] = numpy.copy(counts)
        fc = hifive.FiveC(fivec_fnames["%s_Nora" % meth])
        temp = fc.cis_heatmap(0, datatype='fragment', arraytype='compact', binsize=0, skipfiltered=True)
        ratio2 = temp.shape[1] / float(temp.shape[0])
        imgs["%s_Nora" % meth] = hifive.plotting.plot_full_array(temp, symmetricscaling=False)
        counts = temp[:, :, 0].ravel()
        expected = temp[:, :, 1].ravel()
        if meth == 'Prob':
            temp1 = numpy.zeros((temp.shape[0], temp.shape[1]), dtype=numpy.float32)
            temp1[numpy.where(temp[:, :, 0] > 0.0)] = 1
            imgs["Raw_Nora"] = hifive.plotting.plot_full_array(
                            numpy.dstack((temp[:, :, 0], temp1)), symmetricscaling=False)
            data["Raw_Nora"] = numpy.copy(counts)
            fragments = fc.frags['fragments'][...]
            regions = fc.frags['regions'][...]
            binbounds = numpy.hstack((
                    fragments['start'][regions['start_frag'][0]:regions['stop_frag'][0]].reshape(-1, 1),
                    fragments['stop'][regions['start_frag'][0]:regions['stop_frag'][0]].reshape(-1, 1)))
            binbounds = binbounds[numpy.where(fc.filter[regions['start_frag'][0]:regions['stop_frag'][0]])[0], :]
            temp = hm_nora[0]
            strands = fragments['strand'][regions['start_frag'][0]:regions['stop_frag'][0]]
            temp = temp[numpy.where(strands==0)[0], :, :][:, numpy.where(strands == 1)[0], :]
            imgs["HiC_Nora"] = hifive.plotting.plot_full_array(temp, symmetricscaling=False)
            hic_counts = temp[:, :, 0].ravel()
            hic_expected = temp[:, :, 1].ravel()
            where = numpy.where(hic_expected > 0.0)[0]
            hic_counts[where] /= hic_expected[where]
            data["HiC_Nora"] = numpy.copy(hic_counts)
        where = numpy.where(expected > 0.0)[0]
        counts[where] /= expected[where] 
        data["%s_Nora" % meth] = numpy.copy(counts)
    correlations = {}
    output = open(out_fname.replace('pdf', 'txt'), 'w')
    print >> output, "Method\tPhillips\tNora"
    for meth in methods:
        temp = [meth]
        for name in ["Phillips", "Nora"]:
            valid = numpy.where((data["%s_%s" % (meth, name)] > 0.0) * (data["HiC_%s" % name] > 0.0))
            correlations["%s_%s" % (meth, name)] = numpy.corrcoef(numpy.log(data["%s_%s" % (meth, name)][valid]),
                                                              numpy.log(data["HiC_%s" % name][valid]))[0, 1]
            temp.append(str(correlations["%s_%s" % (meth, name)]))
        print >> output, '\t'.join(temp)
    output.close()
    width = 16.8
    spacer = 0.3
    c = canvas.canvas()
    plot_width = (width - spacer * 3.0 - 0.4) / 4.0
    for i, meth in enumerate(["Raw", "Prob", "HiC"]):
        meth_names = {"Raw":"Raw", "Prob":"HiFive", "HiC":"HiC"}
        c.text(plot_width * (i + 1.5) + spacer * (i + 1), (ratio1 + ratio2) * plot_width + spacer + 0.1,
               "%s" % meth_names[meth], [text.halign.center, text.valign.bottom, text.size(-2)])
        c.insert(bitmap.bitmap(0, 0, imgs["%s_Phillips" % meth], width=plot_width),
                 [trafo.translate((i + 1) * (plot_width + spacer), plot_width * ratio2 + spacer)])
        c.insert(bitmap.bitmap(0, 0, imgs["%s_Nora" % meth], width=plot_width),
                 [trafo.translate((i + 1) * (plot_width + spacer), 0)])
    g = graph.graphxy(width=plot_width - 0.8, height=plot_width * ratio1,
                      x=graph.axis.nestedbar(painter=graph.axis.painter.bar(nameattrs=None)),
                      y=graph.axis.lin(painter=painter),
                      x2=graph.axis.lin(parter=None, min=0, max=1),
                      y2=graph.axis.lin(parter=None, min=0, max=1))
    for i, meth in enumerate(methods):
        Y = numpy.zeros(2, dtype=numpy.float32)
        col = method_colors[meth]
        for j, name in enumerate(["Phillips", "Nora"]):
            Y[j] = correlations["%s_%s" % (meth, name)]
        g.plot(graph.data.points(zip(zip(range(Y.shape[0]), [i] * Y.shape[0]), Y), xname=1, y=2),
               [graph.style.changebar([col])])
    g.text(-0.8, plot_width * ratio1 * 0.5, "Correlation",
           [text.halign.center, text.valign.top, text.size(-3), trafo.rotate(90)])
    g.text((plot_width - 0.8) * 0.25, -0.1, "Phillips",
           [text.halign.center, text.valign.top, text.size(-3)])
    g.text((plot_width - 0.8) * 0.75, -0.1, "Nora",
           [text.halign.center, text.valign.top, text.size(-3)])
    c.insert(g, [trafo.translate(0.8, plot_width * ratio2 + spacer)])
    c.text(width, (ratio1 + ratio2 * 0.5) * plot_width + spacer, "Phillips",
           [text.halign.center, text.valign.top, trafo.rotate(-90), text.size(-2)])
    c.text(width, ratio1 * 0.5 * plot_width, "Nora",
           [text.halign.center, text.valign.top, trafo.rotate(-90), text.size(-2)])
    meth_names = {"Raw":"Raw", "Prob":"HiFive-Probability", "Exp":"HiFive-Express", "Bin":"HiFive-Binning",
                  "Exp-KR":"HiFive-ExpressKR", "Exp-KR-dist":"HiFive-ExpressKR-dist"}
    for i, meth in enumerate(methods):
        c.fill(path.rect(1.0, plot_width * ratio1 - 1.0 - i * 0.5, 0.2, 0.2), [method_colors[meth]])
        c.text(1.3, plot_width * ratio1 - 0.9 - i * 0.5, "%s" % meth_names[meth],
               [text.halign.left, text.valign.middle, text.size(-3)])
    c.writePDFfile(out_fname)


def dynamically_bin(hic1, hic2, region):
    counts = hic1['%i.counts' % region][...] + hic2['%i.counts' % region][...]
    expected = hic1['%i.expected' % region][...] + hic2['%i.expected' % region][...]
    upper = numpy.zeros((counts.shape[0], 2), dtype=numpy.float32)
    upper[:, 0] = counts
    upper[:, 1] = expected
    mids1 = hic1['%i.mids' % region][...]
    mids2 = hic2['%i.mids' % region][...]
    indices = numpy.triu_indices(mids1.shape[0], 1)
    unbinned_counts1 = numpy.zeros((mids1.shape[0], mids1.shape[0]), dtype=numpy.int32)
    unbinned_counts1[indices] = hic1['%i.unbinned_counts' % region][...]
    unbinned_counts1[indices[1], indices[0]] = unbinned_counts1[indices]
    unbinned_expected1 = numpy.zeros((mids1.shape[0], mids1.shape[0]), dtype=numpy.int32)
    unbinned_expected1[indices] = hic1['%i.unbinned_expected' % region][...]
    unbinned_expected1[indices[1], indices[0]] = unbinned_expected1[indices]
    indices = numpy.triu_indices(mids2.shape[0], 1)
    unbinned_counts2 = numpy.zeros((mids2.shape[0], mids2.shape[0]), dtype=numpy.int32)
    unbinned_counts2[indices] = hic2['%i.unbinned_counts' % region][...]
    unbinned_counts2[indices[1], indices[0]] = unbinned_counts2[indices]
    unbinned_expected2 = numpy.zeros((mids2.shape[0], mids2.shape[0]), dtype=numpy.int32)
    unbinned_expected2[indices] = hic2['%i.unbinned_expected' % region][...]
    unbinned_expected2[indices[1], indices[0]] = unbinned_expected2[indices]
    bounds = hic1['%i.bounds' % region][...]
    allmids = numpy.zeros((mids1.shape[0] + mids2.shape[0], 2), dtype=numpy.int32)
    allmids[:mids1.shape[0], 0] = mids1 - 1
    allmids[:mids1.shape[0], 1] = mids1 + 1
    allmids[mids1.shape[0]:, 0] = mids2 - 1
    allmids[mids1.shape[0]:, 1] = mids2 + 1
    allmids = allmids[numpy.argsort(allmids[:, 0]), :]
    indices1 = numpy.searchsorted(allmids[:, 1], mids1)
    indices1_1 = (indices1.reshape(-1, 1) * allmids.shape[0] + indices1.reshape(1, -1)).ravel()
    indices2 = numpy.searchsorted(allmids[:, 1], mids2)
    indices2_1 = (indices2.reshape(-1, 1) * allmids.shape[0] + indices2.reshape(1, -1)).ravel()
    unbinned = numpy.zeros((allmids.shape[0], allmids.shape[0], 2), dtype=numpy.float32)
    unbinned[:, :, 0] += numpy.bincount(indices1_1, minlength=allmids.shape[0] ** 2,
                                        weights=unbinned_counts1.ravel()).reshape(allmids.shape[0], -1)
    unbinned[:, :, 1] += numpy.bincount(indices1_1, minlength=allmids.shape[0] ** 2,
                                        weights=unbinned_expected1.ravel()).reshape(allmids.shape[0], -1)
    unbinned[:, :, 0] += numpy.bincount(indices2_1, minlength=allmids.shape[0] ** 2,
                                        weights=unbinned_counts2.ravel()).reshape(allmids.shape[0], -1)
    unbinned[:, :, 1] += numpy.bincount(indices2_1, minlength=allmids.shape[0] ** 2,
                                        weights=unbinned_expected2.ravel()).reshape(allmids.shape[0], -1)
    indices = numpy.triu_indices(unbinned.shape[0], 1)
    unbinned = unbinned[indices[0], indices[1], :]
    indices = numpy.triu_indices(bounds.shape[0], 1)
    hifive.hic_binning.dynamically_bin_cis_array(unbinned, allmids, upper, bounds,
                                                 expansion_binsize=0, minobservations=25)
    binned = numpy.zeros((bounds.shape[0], bounds.shape[0], 2), dtype=numpy.float32)
    binned[indices[0], indices[1], :] = upper
    binned[indices[1], indices[0], :] = upper
    return binned

if __name__ == "__main__":
    main()
