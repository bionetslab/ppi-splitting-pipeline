"""intact.zip downloaded from https://www.ebi.ac.uk/intact/download/ftp"""
import pandas as pd

intact = pd.read_csv("~/Downloads/intact/intact.txt", sep = "\t")
intact = intact[["#ID(s) interactor A", "ID(s) interactor B", "Interaction detection method(s)", "Taxid interactor A", "Taxid interactor B", "Interaction type(s)"]]
intact_filtered = intact[intact["#ID(s) interactor A"].str.startswith("uniprotkb")]
intact_filtered = intact_filtered[intact_filtered["ID(s) interactor B"].str.startswith("uniprotkb")]
intact_filtered.columns = ["protein1", "protein2", "method", "taxid1", "taxid2", "type"]

intact_filtered["protein1"] = intact_filtered["protein1"].str.replace("uniprotkb:", "")
intact_filtered["protein2"] = intact_filtered["protein2"].str.replace("uniprotkb:", "")
intact_filtered["method"] = intact_filtered["method"].str.replace("psi-mi:", "")
intact_filtered["type"] = intact_filtered["type"].str.replace("psi-mi:", "")
intact_filtered["taxid1"] = intact_filtered["taxid1"].str.replace("taxid:", "")
intact_filtered["taxid1"] = intact_filtered["taxid1"].str.split("|").str[1]
intact_filtered["taxid2"] = intact_filtered["taxid2"].str.replace("taxid:", "")
intact_filtered["taxid2"] = intact_filtered["taxid2"].str.split("|").str[1]

intact_filtered.to_csv("../../data/intact.csv", index=False)
