from dash import Dash, html, dcc, Input, Output, State
import dash_bootstrap_components as dbc
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import base64
from io import BytesIO
from dash.exceptions import PreventUpdate
import os
from gtfparse import read_gtf

gtf = read_gtf(f'Astatotilapia_calliptera.fAstCal1.2.110.gtf')

#top 10 highest peaks included as presets
top_peaks_list = [('chr1', (5919136, 5944136)), ('chr1', (39450585, 40418855)),
                  ('chr3', (3327963, 3352963)), ('chr6', (6492237, 6572923)),
                  ('chr7', (10309548, 10334548)), ('chr7', (44200000, 44225000)),
                  ('chr8', (4074699, 4099699)), ('chr12', (4859506, 4884506)),
                  ('chr12', (5525760, 5550760)), ('chr12', (6053269, 6078269)),
                  ('chr20', (15026647, 15051647)), ('chr23', (20811093, 20836093))]

app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])


app.layout = dbc.Container([
    html.H1("SNP peak explorer", className='mb-2', style={'textAlign': 'center'}),
    dbc.Row([
        dbc.Col([
            dcc.Dropdown(
                id='preset_peaks',
                value='',
                clearable=True,
                options=[{'label': f'{peak[0]}:{peak[1][0]}-{peak[1][1]}', 'value': f'{peak[0]}:{peak[1][0]}-{peak[1][1]}'}
                         for peak in top_peaks_list]
            )
        ], width=4),
        dbc.Col([
            dcc.Dropdown(
                id='chromosome',
                value='chr1',
                clearable=False,
                options=[{'label': f'chr{i}', 'value': f'chr{i}'} for i in range(1, 24) if i != 21]
            )
        ], width=1),
        *[
            dbc.Col([
                dcc.Input(
                    placeholder=f'region {"start" if idx == 0 else "end"}',
                    id=f'region_{"start" if idx == 0 else "end"}',
                    type='number',
                    value=0 if idx == 0 else 25000,
                    step=25000,
                    style={'outline-color': 'yellow'},
                ),
            ], width=2) for idx in range(2)
        ],
        dbc.Col([
            html.A(html.Button('View', id='view-button', className='btn btn-primary'), style={'padding': 20, 'marginTop': 5, 'marginLeft': 2}),
            html.A(html.Button('Download PNG', id='download-png-button', className='btn btn-primary'), style={'padding': 20, 'marginTop': 5, 'marginLeft': 2}, id='download-png-link', href="", download="plot.png")
        ], width=3),
    ]),
    html.Div(id='output-message'),
    dbc.Row([
        dbc.Col([
            html.Img(id='graph-matplotlib')
        ], width=14)
    ])
])

#update inputs when preset peak is selected
@app.callback(
    [Output('chromosome', 'value'), Output('region_start', 'value'), Output('region_end', 'value'), Output('chromosome', 'options')],
    [Input('preset_peaks', 'value')]
)
def update_inputs(selected_peak):
    if selected_peak:
        chromosome, region = selected_peak.split(':')
        region_start, region_end = map(int, region.split('-'))
        chromosome_options = [{'label': f'chr{i}', 'value': f'chr{i}'} for i in range(1, 24) if i != 21]
        return chromosome, region_start, region_end, chromosome_options
    return 'chr1', 0, 25000, [{'label': f'chr{i}', 'value': f'chr{i}'} for i in range(1, 24) if i != 21]

#download PNG of matplotlib figure as produced
@app.callback(
    Output("download-png-link", "href"),
    [Input("download-png-button", "n_clicks")],
    [State('chromosome', 'value'), State('region_start', 'value'), State('region_end', 'value')]
)
def download_png(n_clicks, chromosome, region_start, region_end):
    if not n_clicks:
        raise PreventUpdate
    file_chrom_path = f"data/out_{chromosome}.tsv"
    fig_src, messages = triple_plots(chromosome, region_start, region_end, gtf, file_chrom_path)
    plot_data = fig_src.split(",")[1]
    png_filename = f"plot_{chromosome}_{region_start}_{region_end}.png"
    png_filepath = os.path.join(png_filename)
    with open(png_filepath, "wb") as f:
        f.write(base64.b64decode(plot_data))
    return png_filepath

#Update plot, when view button selected with selected regions
@app.callback(
    [Output('graph-matplotlib', 'src'), Output('output-message', 'children')],
    [Input('view-button', 'n_clicks')],
    [State('chromosome', 'value'), State('region_start', 'value'), State('region_end', 'value')]
)
def update_plot(n_clicks, chromosome, region_start, region_end):
    messages = []
    file_chrom_path = f"data/out_{chromosome}.tsv"
    if n_clicks > 0:
        fig_src, additional_messages = triple_plots(chromosome, region_start, region_end, gtf, file_chrom_path)
        messages.extend(additional_messages)
        return fig_src, messages
    return None, messages

#plotting function
def triple_plots(chromosome, region_start, region_end, gtf, file_chrom_path):
    messages = []
    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(14, 12), gridspec_kw={'height_ratios': [1.5, 1, 1]}, sharex=True)
    axGENE, axHSCAN, axPBS = axs
    gtf = gtf[((gtf['feature'] == 'exon') | (gtf['feature'] == 'gene')) & (gtf['seqname'] == chromosome.replace("chr", ""))]
    merged_df = pd.read_csv(file_chrom_path, sep='\t', index_col=[0])
    gene_starts = gtf[gtf['feature'] == 'gene']['start']
    gene_ends = gtf[gtf['feature'] == 'gene']['end']
    overlapping_genes = ((gene_starts <= region_start) & (gene_ends >= region_start)) | ((gene_starts <= region_end) & (gene_ends >= region_end))
    if overlapping_genes.any():
        for item in [("start", region_start), ("end", region_end)]:
            messages.append(html.Div(f"The specified region overlaps with gene(s) on {chromosome}. "
                                      f"The region {item[0]} has been shifted from {item[1]} to "
                                      f"{min(item[1], gene_starts[overlapping_genes].min())}", style={'color': 'red'}))
        region_start = int(min(region_start, gene_starts[overlapping_genes].min()))
        region_end = int(max(region_end, gene_ends[overlapping_genes].max()))
    gtf = gtf[(gtf['start'] >= region_start) & (gtf['end'] <= region_end)]
    gene_groups = gtf.groupby(gtf['gene_id'])
    num_genes = len(gene_groups)
    ylim = (-0.1, (num_genes + 1) / num_genes) if num_genes > 0 else (0, 1)
    axGENE.set_ylim(ylim)
    for ax in axs:
        ax.set_xlim(region_start, region_end)
    arrow_head_size = (region_end - region_start) / 150
    for i, (name, group) in enumerate(gene_groups):
        gene_name = str(group['gene_name'].iloc[0])
        label = gene_name if gene_name != "" else name
        strand = group.iloc[0]['strand']
        gene_start = group[group['feature'] == 'gene'].iloc[0]['start']
        gene_end = group[group['feature'] == 'gene'].iloc[0]['end']
        exon_starts = group[group['feature'] == 'exon']['start'].tolist()
        exon_ends = group[group['feature'] == 'exon']['end'].tolist()
        exon_starts.sort()
        exon_ends.sort()
        y_pos = (i + 0.2) / num_genes
        for j in range(len(exon_starts)):
            if j > 0:
                axGENE.plot([exon_ends[j - 1], exon_starts[j]], [y_pos, y_pos], color='black', linewidth=1)
            axGENE.fill_between([exon_starts[j], exon_ends[j]], [y_pos + 0.05, y_pos + 0.05], [y_pos - 0.05, y_pos - 0.05], color='gray')
        if strand == '+':
            arrow_length = gene_end - gene_start
            arrow_end = gene_end - arrow_head_size if arrow_length > 2 * arrow_head_size else gene_end
            axGENE.arrow(gene_start, y_pos + 0.075, arrow_end - gene_start, 0, head_width=0.05, head_length=arrow_head_size, fc='k', ec='k', label='Gene')
        elif strand == '-':
            arrow_length = gene_end - gene_start
            arrow_end = gene_start + arrow_head_size if arrow_length > 2 * arrow_head_size else gene_start
            axGENE.arrow(gene_end, y_pos + 0.075, arrow_end - gene_end, 0, head_width=0.05, head_length=arrow_head_size, fc='k', ec='k', label='Gene')
        if num_genes <= 4:
            axGENE.text((gene_start + gene_end) / 2, y_pos + 0.13, label, ha='center', va='center')
        if 4 < num_genes < 20:
            axGENE.text((gene_start + gene_end) / 2, y_pos + 0.13, label, ha='right', va='center')
    axGENE.set_ylim(ylim)
    axGENE.set_yticks([])
    axGENE.set_ylabel("Ensembl Genes", fontsize=14)
    # axGENE.spines['left'].set_visible(False)
    # axGENE.spines['right'].set_visible(False)
    # axGENE.spines['top'].set_visible(False)
    for ax, col, ylabel in zip([axHSCAN, axPBS], ['H', 'PBS_pop1'], ['H value', 'PBS score']):
        merged_df[col].plot(x='pos', y=col, style='.', rasterized=(ax == axHSCAN), color='#3BC14A' if ax == axHSCAN else 'blue', ax=ax)
        ax.set_ylabel(ylabel, fontsize=14)
        ax.set_xlabel('position', fontsize=14)
    buf = BytesIO()
    fig.savefig(buf, format="png")
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,{fig_data}'
    return fig_bar_matplotlib, messages

if __name__ == '__main__':
    app.run_server(debug=False, port=8002)
