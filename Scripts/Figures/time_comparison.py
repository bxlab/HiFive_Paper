#!/usr/bin/env python

import sys
import os
from glob import glob
from math import floor

import numpy
from pyx import canvas, text, path, graph, color, trafo, unit, attr, deco, style, bitmap


unit.set(defaultunit="cm")
text.set(mode="latex")
text.preamble(r"\usepackage{times}")
text.preamble(r"\usepackage{sansmath}")
text.preamble(r"\sansmath")
text.preamble(r"\renewcommand*\familydefault{\sfdefault}")
painter = graph.axis.painter.regular( labeldist=0.1, labelattrs=[text.size(-3)], titleattrs=[text.size(-3)] )
painter2 = graph.axis.painter.regular( labelattrs=None )

step_colors = {
    0:color.cmyk.MidnightBlue,
    1:color.cmyk.Mahogany,
    2:color.cmyk.CadetBlue,
    3:color.cmyk.Dandelion,
    4:color.cmyk.OliveGreen,
}

step_names = {
    '0':'Loading data',
    '1':'Filtering and preprocessing',
    '2':'Binning',
    '3':'Normalization',
    '4':'Heatmapping'
}

def main():
    width = 16.8
    out_fname = sys.argv[1]
    basedir = "%s/Analysis/Timing" % '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-2])
    data_fnames = {
        "HiFive-Probability":{'0':"%s/hifive_data" % basedir, '1':"%s/hifive_project" % basedir,
                '3':'%s/hifive_prob' % basedir, '4':"%s/hifive_prob_heatmap" % basedir},
        "HiFive-Binning":{'0':"%s/hifive_data" % basedir, '1':"%s/hifive_project_nodist" % basedir,
               '3':'%s/hifive_bin' % basedir, '4':"%s/hifive_bin_heatmap" % basedir},
        "HiFive-Express":{'0':"%s/hifive_data" % basedir, '1':"%s/hifive_project" % basedir,
               '3':'%s/hifive_exp' % basedir, '4':"%s/hifive_exp_heatmap" % basedir},
        "HiFive-ExpressKR":{'0':"%s/hifive_data" % basedir, '1':"%s/hifive_project_nodist" % basedir,
                 '3':'%s/hifive_expKR' % basedir, '4':"%s/hifive_expKR_heatmap" % basedir},
        "HiFive-ExpressKR w/distance":{'0':"%s/hifive_data" % basedir, '1':"%s/hifive_project" % basedir,
                     '3':'%s/hifive_expKRdist' % basedir, '4':"%s/hifive_expKRdist_heatmap" % basedir},
        "HiCPipe":{'0':"%s/bam2raw" % basedir, '1':"%s/hicpipe_data" % basedir,
                   '2':'%s/hicpipe_binning' % basedir, '3':'%s/hicpipe_norm' % basedir,
                   '4':'%s/hicpipe_heatmap' % basedir},
        "HiCLib":{'0':"%s/hiclib_mapping" % basedir, '1':'%s/hiclib_data' % basedir,
                  '3':"%s/hiclib_norm" % basedir, '4':'%s/hiclib_heatmap' % basedir},
        "HiCNorm":{'0':"%s/bam2raw" % basedir, '1':"%s/hicpipe_data" % basedir,
                   "2":"%s/hicnorm_data" % basedir, '3':'%s/hicnorm_norm' % basedir},
    }
    data = load_data(data_fnames)
    c = canvas.canvas()
    c.insert(plot_bargraph(data, width, 4.0), [trafo.translate(4.0, 0)])
    c.insert(plot_key(width * 0.3, 1.5), [trafo.translate(width * 0.75 - 1.0, 0.2)])
    c.writePDFfile(out_fname)

def plot_key(width, height):
    c = canvas.canvas()
    wstep = width / 2.0
    hstep = height / 3
    for i in range(5):
        x = wstep * (i / 3)
        y = height - (i % 3 + 0.5) * hstep
        c.fill(path.rect(x + 0.1, y - 0.15, 0.3, 0.3), [step_colors[i]])
        temp = step_names[str(i)].split(' ')
        if len(temp) > 2:
            c.text(x + 0.5, y + 0.13, ' '.join(temp[:2]),
                   [text.halign.left, text.valign.middle, text.size(-2)])
            c.text(x + 0.5, y - 0.13, temp[-1],
                   [text.halign.left, text.valign.middle, text.size(-2)])
        else:
            c.text(x + 0.5, y, step_names[str(i)],
                   [text.halign.left, text.valign.middle, text.size(-2)])
    return c

def plot_bargraph(data, width, height):
    methods = ['HiCLib', 'HiCPipe', 'HiCNorm', 'HiFive-Probability', 'HiFive-Binning', 'HiFive-Express',
               'HiFive-ExpressKR', 'HiFive-ExpressKR w/distance']
    ho = 4.0
    left_width = (width - ho) * 0.45
    mid_width1 = (width - ho) * 0.3
    mid_width2 = (width - ho) * 0.125
    right_width = (width - ho) * 0.125
    bar_height = height / len(methods) - 0.1
    data_totals = {}
    ranges = numpy.zeros((4, 2), dtype=numpy.float32)
    for meth in data:
        data_totals[meth] = find_total(data[meth])
        if meth == 'HiCPipe':
            ranges[1, 1] = data_totals[meth]
        elif meth == 'HiCNorm':
            ranges[2, 1] = data_totals[meth]
        elif meth == 'HiFive-Probability':
            ranges[3, 1] = data_totals[meth]
        else:
            ranges[0, 1] = max(ranges[0, 1], data_totals[meth])
    ranges /= 60.0
    ranges[0, 1] = 28.0
    ranges[1, 0] = ranges[1, 1] - ranges[0, 1] / 0.45 * 0.3 * 0.9
    ranges[1, 1] = ranges[1, 1] + ranges[0, 1] / 0.45 * 0.3 * 0.1
    ranges[2, 0] = ranges[2, 1] - ranges[0, 1] / 0.45 * 0.125 * 0.5
    ranges[2, 1] = ranges[2, 1] + ranges[0, 1] / 0.45 * 0.125 * 0.5
    ranges[3, 0] = ranges[3, 1] - ranges[0, 1] / 0.45 * 0.125 * 0.5
    ranges[3, 1] = ranges[3, 1] + ranges[0, 1] / 0.45 * 0.125 * 0.5
    c = canvas.canvas()
    g1 = graph.graphxy(width=left_width, height=height,
                       x=graph.axis.lin(painter=painter, min=0, max=ranges[0, 1]),
                       x2=graph.axis.lin(parter=None, min=0, max=ranges[0, 1]),
                       y=graph.axis.lin(parter=None, min=0, max=1),
                       y2=graph.axis.lin(painter=None, min=0, max=1))
    c.insert(g1)
    g2 = graph.graphxy(width=mid_width1, height=height,
                       x=graph.axis.lin(painter=painter, min=ranges[1, 0], max=ranges[1, 1]),
                       x2=graph.axis.lin(parter=None, min=ranges[1, 0], max=ranges[1, 1]),
                       y2=graph.axis.lin(painter=None, min=0, max=1),
                       y=graph.axis.lin(painter=None, min=0, max=1))
    c.insert(g2, [trafo.translate(left_width, 0)])
    g3 = graph.graphxy(width=mid_width2, height=height,
                       x=graph.axis.lin(painter=painter, min=ranges[2, 0], max=ranges[2, 1]),
                       x2=graph.axis.lin(parter=None, min=ranges[2, 0], max=ranges[2, 1]),
                       y2=graph.axis.lin(painter=None, min=0, max=1),
                       y=graph.axis.lin(painter=None, min=0, max=1))
    c.insert(g3, [trafo.translate(left_width + mid_width1, 0)])
    g4 = graph.graphxy(width=right_width, height=height,
                       x=graph.axis.lin(painter=painter, min=ranges[3, 0], max=ranges[3, 1]),
                       x2=graph.axis.lin(parter=None, min=ranges[3, 0], max=ranges[3, 1]),
                       y2=graph.axis.lin(parter=None, min=0, max=1),
                       y=graph.axis.lin(painter=None, min=0, max=1))
    c.insert(g4, [trafo.translate(left_width + mid_width1 + mid_width2, 0)])
    split = canvas.canvas()
    split.fill(path.path(path.moveto(-0.15, -0.2), path.lineto(0.05, 0.2), path.lineto(.15, 0.2),
               path.lineto(-0.05, -0.2), path.closepath()), [color.cmyk.White])
    split.stroke(path.line(-0.15, -0.2, 0.05, 0.2))
    split.stroke(path.line(-0.05, -0.2, 0.15, 0.2))
    c.insert(split, [trafo.translate(left_width, 0)])
    c.insert(split, [trafo.translate(left_width, height)])
    c.insert(split, [trafo.translate(left_width + mid_width1, 0)])
    c.insert(split, [trafo.translate(left_width + mid_width1, height)])
    c.insert(split, [trafo.translate(left_width + mid_width1 + mid_width2, 0)])
    c.insert(split, [trafo.translate(left_width + mid_width1 + mid_width2, height)])
    for i, meth in enumerate(methods):
        c.insert(plot_bar(data[meth], ranges, bar_height, left_width / ranges[0, 1], split),
                 [trafo.translate(0, height - 0.05 - bar_height * (i + 1) - i * 0.1)])
        c.text(-0.1, height * (len(methods) - i - 0.5) / len(methods), meth,
               [text.halign.right, text.valign.middle, text.size(-2)])
    c.text((width - ho) / 2.0, -0.35, "Runtime (minutes)",
           [text.halign.center, text.valign.top, text.size(-2)])
    return c

def plot_bar(data, ranges, height, scale, split):
    c = canvas.canvas()
    pos = 0.0
    x1 = 0.0
    for i in range(5):
        if str(i) in data:
            x2 = data[str(i)] / 60.0 + pos
            pos2 = x2
            for j in range(1, ranges.shape[0]):
                if x2 > ranges[j, 0] and x2 < ranges[j, 1]:
                    for k in range(j):
                        x2 -= ranges[k + 1, 0] - ranges[k, 1]
                    break
            c.fill(path.rect(x1 * scale, 0, (x2 - x1) * scale, height), [step_colors[i]])
            split_pos = ranges[0, 1] * scale
            for j in range(ranges.shape[0] - 1):
                if pos2 > ranges[j, 1]:
                    c.insert(split, [trafo.translate(split_pos, height * 0.5)])
                    split_pos += (ranges[j + 1, 1] - ranges[j + 1, 0]) * scale
            pos += data[str(i)] / 60.0
            x1 = x2
    return c

def find_max(data):
    max_time = 0.0
    for name in data:
        max_time = max(max_time, find_total(data[name]))
    return max_time

def find_total(data):
    time = 0.0
    for key in data:
        time += data[key]
    return time

def load_data(methods):
    data = {}
    for name in methods:
        data[name] = {}
        for key in methods[name]:
            data[name][key] = find_median(methods[name][key], 'Elapsed')
    return data

def find_median(fprefix, name):
    fnames = glob("%s_[0-9].txt" % fprefix)
    if len(fnames) == 0:
        return 0.0
    times = []
    for fname in fnames:
        for line in open(fname):
            temp = line[1:-1].split()
            if len(temp) == 0:
                continue
            if temp[0] == name and name == 'Elapsed':
                temp1 = temp[-1].split(':')
                time = float(temp1[-1])
                for i in range(len(temp1) - 1)[::-1]:
                    time += 60.0 ** (len(temp1) - i - 1) * int(temp1[i])
                times.append(time)
            elif temp[0] == name:
                times.append(float(temp[-1]) / 1000.0)
    if len(times) > 0:
        times = numpy.array(times)
        return numpy.median(times)
    else:
        return 0.0

if __name__ == "__main__":
    main()
