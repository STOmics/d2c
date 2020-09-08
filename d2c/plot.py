# -*- coding: utf-8 -*-

import os
import sys
import math
import gzip
import time

from numpy import log10
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import plotly.figure_factory as ff
from plotly.subplots import make_subplots

DEBUG = False


def scatterBarcode(filename, threshold, out_prefix):
    
    # barcode_counts = []
    # with open(filename) as fh_in:
    #     for line in fh_in:
    #         line = line.strip().split(',')
    #         barcode_counts.append(int(line[1]))

    # barcode_counts.sort(reverse=True)
    start_time = time.time()
    df = pd.read_csv(filename, names = ["barcode", "cnt"])
    df.sort_values(by="cnt", ascending=False, ignore_index=True)
    
    # Plot bead barcode knee
    # df = DataFrame(barcode_counts, index=list(range(0,len(barcode_counts))), columns=["cnt"])
    df['PassKnee'] = df["cnt"] > threshold

    if DEBUG:
        end_time = time.time()
        print(f"repare data time(s): {end_time-start_time}")
        start_time = time.time()

    fig = px.scatter(
        df,
        x=df.index,
        y="cnt",
        log_x=True,
        log_y=True,
        color="PassKnee",
        labels={"index": "Barcode in rank-descending order", "cnt":"Reads per barcode"}
    )
    
    fig.update_layout(title_text="BeadBarcodeKnee", title_x=0.5)
    
    x_intercept_bead = len(df['cnt'][df['cnt'] > threshold])
    fig.add_shape(
        type="line", line_color="LightSeaGreen", line_width=3, opacity=1,
        x0=0, x1=1, xref="paper", y0=threshold, y1=threshold, yref="y"
    )
    fig.add_shape(
        type="line", line_color="LightSeaGreen", line_width=3, opacity=1,
        x0=x_intercept_bead, x1=x_intercept_bead, xref="x", y0=0, y1=1, yref="paper"
    )
    if DEBUG:
        end_time = time.time()
        print(f"plot barcode knee time(s): {end_time-start_time}")
        start_time = time.time()

    fig.write_image(out_prefix+".BeadBarcodeKnee.png")
    if DEBUG:
        end_time = time.time()
        print(f"dump barcode knee time(s): {end_time-start_time}")
        start_time = time.time()
    fig.write_html(out_prefix+".BeadBarcodeKnee.html")
    if DEBUG:
        end_time = time.time()
        print(f"dump barcode knee time(s): {end_time-start_time}")
        start_time = time.time()

    # Plot bead barcode knee density
    df2 = log10(df[df["cnt"] > 50]["cnt"]+1)
    fig = ff.create_distplot([df2], ['density'], show_hist=False, show_rug=False,
                        curve_type='kde')
                        
    fig.update_xaxes(title="Count per barcode log10")
    fig.update_yaxes(title="Density")
    fig.update_layout(title_text="BeadBarcodeKneeDensity", title_x=0.5, showlegend=False)
    
    x_intercept_bead_log = math.log10(x_intercept_bead)
    fig.add_shape(
            type="line", line_width=3, opacity=1,
            x0=x_intercept_bead_log, x1=x_intercept_bead_log, xref="x", y0=0, y1=1, yref="paper"
        )

    fig.write_image(out_prefix+".BeadBarcodeKneeDensity.png")
    fig.write_html(out_prefix+".BeadBarcodeKneeDensity.html")
    if DEBUG:
        end_time = time.time()
        print(f"plot density time(s): {end_time-start_time}")
        start_time = time.time()

    # Plot bead barcode knee curve
    df["cumsum"] = df["cnt"].cumsum()
    fig = px.line(df, x=df.index, y="cumsum")
    fig.add_shape(
        type="line", line_color="LightSeaGreen", line_width=1, opacity=1,
        x0=threshold, x1=threshold, xref="x", y0=0, y1=1, yref="paper"
    )
    fig.update_xaxes(title="Barcode rank")
    fig.update_yaxes(title="Cumulative read count")
    fig.update_layout(title_text="BeadBarcodeKneeCurve", title_x=0.5, showlegend=False)

    if DEBUG:
        end_time = time.time()
        print(f"plot curve time(s): {end_time-start_time}")
        start_time = time.time()

    fig.write_image(out_prefix+".BeadBarcodeKneeCurve.png")
    if DEBUG:
        end_time = time.time()
        print(f"dump curve time(s): {end_time-start_time}")
        start_time = time.time()
    fig.write_html(out_prefix+".BeadBarcodeKneeCurve.html")
    if DEBUG:
        end_time = time.time()
        print(f"dump curve time(s): {end_time-start_time}")
        start_time = time.time()
    
def scatterJaccard(filename, threshold, out_prefix):
    barcode_counts = []
    with gzip.open(filename) as fh_in:
        fh_in.readline()
        i = 1
        for line in fh_in:
            line = line.decode().strip().split(',')
            barcode_counts.append(float(line[-2]))
            i += 1
            if i == 1000000:
                break

    barcode_counts.sort(reverse=True)
    df = pd.DataFrame(barcode_counts, index=list(range(0,len(barcode_counts))), columns=["cnt"])
    df['PassKnee'] = df["cnt"] > threshold
    
    fig = px.scatter(
        df,
        x=df.index,
        y="cnt",
        log_x=True,
        log_y=True,
        color="PassKnee",
        labels={"index": "d2c overlap score in rank-descending order", "cnt":"d2c overlap score per barcode pair"}
    )
    fig.update_layout(title_text="JaccardOverlapKnee", title_x=0.5)
    
    x_intercept_bead = len(df['cnt'][df['cnt'] > threshold])
    fig.add_shape(
        type="line", line_color="LightSeaGreen", line_width=3, opacity=1,
        x0=0, x1=1, xref="paper", y0=threshold, y1=threshold, yref="y"
    )
    fig.add_shape(
        type="line", line_color="LightSeaGreen", line_width=3, opacity=1,
        x0=x_intercept_bead, x1=x_intercept_bead, xref="x", y0=0, y1=1, yref="paper"
    )
    
    fig.write_image(out_prefix+".JaccardOverlapKnee.png")
    fig.write_html(out_prefix+".JaccardOverlapKnee.html")
    
    

def scatterSaturation(filename, out_prefix):
    x = [0]
    y = [[0],[0]]
    with open(filename) as fh_in:
        for line in fh_in:
            if line.startswith('#'): continue
            line = line.strip().split()
            x.append(line[1])
            y[0].append(line[2])
            y[1].append(line[3])
            
    fig = make_subplots(rows=1, cols=2,
            subplot_titles=("SequencingSaturation", "MedianUniqueFragmentsPerCell"))
    fig.add_trace(
        go.Scatter(x=x, y=y[0], name="SeqSat"),
        row=1, col=1
    )
    fig.add_trace(
        go.Scatter(x=x, y=y[1], name="Fragments"),
        row=1, col=2
    )

    # Update xaxis properties
    xaxes= "Mean Reads per Cell"
    fig.update_xaxes(title_text=xaxes, row=1, col=1)
    fig.update_xaxes(title_text=xaxes, row=1, col=2)
    
    # Update yaxis properties
    yaxes_saturation = "Sequencing Saturation"
    yaxes_median = "Median Unique Fragments per Cell"
    fig.update_yaxes(title_text=yaxes_saturation, row=1, col=1)
    fig.update_yaxes(title_text=yaxes_median, row=1, col=2)
   
    fig.update_layout(height=800, width=1000, showlegend=False)

    fig.write_image(out_prefix+".SequencingSaturation.png")
    fig.write_html(out_prefix+".SequencingSaturation.html")
    
    
def main(file_path, prefix):
    out_path = os.path.join(file_path, "plot")
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    out_prefix = os.path.join(out_path, prefix)
    runname = os.path.join(file_path, prefix)

    start_time = time.time()
    seq_sat_file = runname+'.sequenceSaturation.tsv'
    print(f'seq_sat_file={seq_sat_file}')
    if os.path.exists(seq_sat_file):
        scatterSaturation(seq_sat_file, out_prefix)
        if DEBUG:
            end_time = time.time()
            print(f"plot saturation time(s): {end_time-start_time}")
            start_time = time.time()

    param_file = runname+'.d2cParam.csv'
    print(f'param_file={param_file}')
    if not os.path.exists(param_file):
        print("Not exists param file: ", param_file)
        sys.exit(1)

    params = {}
    with open(param_file) as fh_in:
        for line in fh_in:
            line = line.strip().split(',')
            params[line[0]] = float(line[1])
            
    bead_file = runname+'.barcodeQuantSimple.csv'
    print(f'bead_file={bead_file}')
    if 'bead_threshold' in params and os.path.exists(bead_file):
        bead_threshold = params['bead_threshold']
        print(f'bead_threshold={bead_threshold}')
        scatterBarcode(bead_file, bead_threshold, out_prefix)
        if DEBUG:
            end_time = time.time()
            print(f"plot barcode time(s): {end_time-start_time}")
            start_time = time.time()
        
    jaccard_file = runname+'.implicatedBarcodes.csv.gz'
    print(f'jaccard_file={jaccard_file}')
    if 'jaccard_threshold' in params and os.path.exists(jaccard_file):
        jaccard_threshold = params['jaccard_threshold']
        print(f'jaccard_threshold={jaccard_threshold}')
        scatterJaccard(jaccard_file, jaccard_threshold, out_prefix)
        if DEBUG:
            end_time = time.time()
            print(f"plot jaccard time(s): {end_time-start_time}")
            start_time = time.time()
        
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Enter <path> <prefix>")
        sys.exit(1)

    file_path = sys.argv[1]
    prefix = sys.argv[2]
    print(f'file_path={file_path} prefix={prefix}')
    main(file_path, prefix)

