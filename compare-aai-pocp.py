import sys
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main(args):
    print('starting POCP vs AAI comparison...')
    # read original pocp table and stack to make long-form
    orig_pocp = 'brady_pocp_table.tab'
    op = pd.read_csv(orig_pocp, sep='\t', index_col=0)
    op_stack = op.stack().reset_index()
    op_stack.columns = ['gA', 'gB', 'POCP']

    # read sourmash AAI table and stack to make long-form
    sourmash_aai = args.sourmash_compare_csv
    sa = pd.read_csv(sourmash_aai, sep=',')#, index_col=0)
    sa.index = sa.columns
    sa_stack = sa.stack().reset_index()
    sa_stack.columns = ['gA', 'gB', 'cAAI']
    #muplitply by 100 to get percentage
    sa_stack['cAAI'] = sa_stack['cAAI'] * 100
    
    #merge sa_stack and op_stack
    merged = pd.merge(sa_stack, op_stack, on=['gA', 'gB'])
    # write an assertion to check that merged has the same number of columns as sa_stack adn op_stack
    assert merged.shape[0] == sa_stack.shape[0] == op_stack.shape[0]
    # remove the diagonal
    merged = merged[merged['gA'] != merged['gB']]
    if args.output_comparison_csv:
        merged.to_csv(args.output_comparison_csv, index=False, sep=',')
        print(f'comparison table saved to {args.output_comparison_csv}')

    # plot scatterplot with seaborn
    sns.set(style='white', rc={'figure.figsize':(11.7,8.27)})
    sns.set_context("paper", font_scale=1.5)
    # scatterplot but dont include labels in legend
    sns.scatterplot(data=merged, x="cAAI", y="POCP", hue="gA", legend=False, palette="Set2")
    # add dashed regression line to plot, print equation above line
    sns.regplot(data=merged, x="cAAI", y="POCP", scatter=False, color="dimgray", 
                line_kws={'linestyle':'--'}, 
                label='y={0:.2f}x+{1:.2f}'.format(*np.polyfit(merged['cAAI'], merged['POCP'], 1)))
    # compute correlation coefficient
    corr = merged['cAAI'].corr(merged['POCP'], method='pearson')
    # print correlation coefficient on bottom right of plot
    plt.text(0.75, 0.85, 'r = {0:.2f}'.format(corr), transform=plt.gca().transAxes, fontsize=14,
                verticalalignment='top', color='dimgray')

    # only show legend for regression line
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., 
            labels=['y={0:.2f}x+{1:.2f}'.format(*np.polyfit(merged['cAAI'], merged['POCP'], 1))])
    plt.savefig(args.output_plot, bbox_inches='tight')
    print(f'plot saved to {args.output_plot}')


# write argparse to input files and parameter
def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument('--pocp-table', default='brady_pocp_table.tab', help='Original POCP table')
    p.add_argument('--sourmash-compare-csv', default='output.pocp/brady.protein.sc5.compare.csv', help='Sourmash compare csv')
    p.add_argument('--output-comparison-csv')
    p.add_argument('--output-plot', default='plots/POCP-vs-cANI.sc5.png', help='Output plot')
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
    
    
    # sourmash_aai = 'output.pocp/brady.protein.sc5.compare'
    # sa = np.load(open(sourmash_aai, 'rb'))
    # labeltext = [x.strip() for x in open(sourmash_aai + '.labels.txt')]
    # print(labeltext[:5])
    # sd = pd.DataFrame(sa, columns = labeltext, index = labeltext)
    # sd.head()
    # load sourmash compare csv (matrix)

