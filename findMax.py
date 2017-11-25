import pandas as pd

def main():
    df = pd.read_csv("colon_data.csv")
    dfDict = df.to_dict(orient = 'index')
    mapdf = pd.read_csv("mapper.csv")
    mapdf.drop('Unnamed: 0', axis=1, inplace = True)
    mDict = mapdf.to_dict('index')
    mapper = dict()

    # make a dictionary that maps each probe to its gene
    for i in mDict.keys():
        entrez_id = mDict[i]['entrez']
        if entrez_id not in mapper.keys():
            mapper[entrez_id] = [mDict[i]['probe']]
        else:
            mapper[entrez_id] +=  [mDict[i]['probe']]
    
    colon_cols = mapper.keys()
    
    new_data = pd.DataFrame(columns= colon_cols)
    
    # for every gene, find the max value amongst its probes
    this_data = []
    for obs in range(df.shape[0]):
        this_row = []
        samps = dfDict[dfDict.keys()[obs]]
        for genes in mapper.keys():
            max_entry = 0
            for i in mapper[genes]:
                col_value = samps[i]
                if col_value > 0:
                    max_entry = col_value
            this_row += [max_entry]
        this_data += [this_row]
    new_data = pd.DataFrame.from_records(this_data, columns = colon_cols)
    
    new_data.to_csv("mod_colon.csv", header = True)
   

main()
