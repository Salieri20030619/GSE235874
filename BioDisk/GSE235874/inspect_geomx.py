import pandas as pd, openpyxl

files = {
    'matrix_q3': 'GSE235874_QC_BQC_Q3Norm_matrixonly.xlsx',
    'seg':       'GSE235874_segmentSummary_readTypes.xlsx',
    'probe':     'GSE235874_BioProbeInformation.xlsx',
}

for k, f in files.items():
    wb = openpyxl.load_workbook(f, read_only=True)
    print(f"\n########## {k}: {f}")
    print("  sheets:", wb.sheetnames)
    for sh in wb.sheetnames:
        ws = wb[sh]
        print(f"    [{sh}] rows={ws.max_row}  cols={ws.max_column}")
    wb.close()

seg = pd.read_excel(files['seg'])
print("\n== segmentSummary 列名 ==")
print(list(seg.columns))
print(seg.head(10).to_string())

mx = pd.read_excel(files['matrix_q3'], nrows=5)
print("\n== Q3 matrix 预览 ==")
print("列数(ROI数≈):", mx.shape[1])
print("前 8 列名:", list(mx.columns[:8]))
print(mx.iloc[:5, :6].to_string())

pb = pd.read_excel(files['probe'], nrows=5)
print("\n== BioProbe 列名 ==")
print(list(pb.columns))
