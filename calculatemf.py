def calculate_mf(dataframe,sample,s_ctx,sum_list = []):
    # for a given ctx, add up all the AD and divide this by the sum of the DP. call this Z_(ctx name)
    # once you have this for every ctx, the P_s for a given ctx is Z_(ctx name)/sum of all Z's
    if 'DP' in dataframe.columns:
        sum_list = []
        dataframe = dataframe[dataframe['SAMPLE']==sample]
        s_ctx_df = dataframe[dataframe['CTX']==s_ctx]
        s_AD = s_ctx_df['ALT_DEPTH'].sum()
        s_DP = s_ctx_df['DP'].sum()
        numerator = s_AD/s_DP
    
        unique_values = dataframe['CTX'].unique().tolist()
        
        for ctx in unique_values:
            s_ctx_df = dataframe[dataframe['CTX']==ctx]
            s_AD = s_ctx_df['ALT_DEPTH'].sum()
            s_DP = s_ctx_df['DP'].sum()
            Z = s_AD/s_DP
            sum_list.append(Z)
        denominator = sum(sum_list)
        return numerator/denominator
    else:
        sum_list = []
        if 'ALT_DEPTH' not in dataframe.columns:
            dff = dataframe['ALT_DEPTH'] = 1
        else:
            dff = dataframe
        dff = dff[dff['SAMPLE']==sample]
        s_ctx_df = dff[dff['CTX']==s_ctx]
        s_AD = s_ctx_df['ALT_DEPTH'].sum()
        numerator = s_AD
    
        unique_values = dataframe['CTX'].unique().tolist()
        
        for ctx in unique_values:
            s_ctx_df = dff[dff['CTX']==ctx]
            s_AD = s_ctx_df['ALT_DEPTH'].sum()
            Z = s_AD
            sum_list.append(Z)
        denominator = sum(sum_list)
        return numerator/denominator
        
# need to make update using groupby as this will be inneficient for large datasets, and Salk might use large datasets...
