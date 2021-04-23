import pandas as pd


def prepare_sleuth_results(sleuth_results):
    df = pd.read_csv(sleuth_results)
    df = df.drop(df.columns[[0,2, 3, 4,7,8,9,10,11,12]], axis=1)
    return df


def prepare_targets(target_file):
    df = pd.read_csv(target_file)
    df = df.drop(df.columns[[1,3,5,7]], axis=1)
    df['6mer'] = df['6mer'].str[:15]
    df['7mera1'] = df['7mera1'].str[:15]
    df['7merm8'] = df['7merm8'].str[:15]
    df['8mer'] = df['8mer'].str[:15]
    return df


def merge(sleuth, targets, output):
    df_temp = targets.drop(targets.columns[[1,2,3]] , axis=1)
    targets_6mer= pd.merge(sleuth, df_temp, left_on = 'target_id', right_on = '6mer', how = 'inner')
    df_temp = targets.drop(targets.columns[[0, 2, 3]], axis=1)
    targets_7mera1= pd.merge(sleuth, df_temp, left_on='target_id', right_on='7mera1', how='inner')

    df_temp = targets.drop(targets.columns[[0, 1, 3]], axis=1)
    targets_7merm8= pd.merge(sleuth, df_temp, left_on='target_id', right_on='7merm8', how='inner')

    df_temp = targets.drop(targets.columns[[0, 1, 2]], axis=1)
    targets_8mer = pd.merge(sleuth, df_temp, left_on='target_id', right_on='8mer', how='inner')

    targets_6mer['6mer'] = '6mer (N= ' + str(len(targets_6mer)) +")"
    targets_6mer.rename(columns={"6mer": "type"}, inplace=True)
    targets_7mera1['7mera1'] = '7mera1 (N= ' + str(len(targets_7mera1)) +")"
    targets_7mera1.rename(columns={"7mera1": "type"}, inplace=True)
    targets_7merm8['7merm8'] = '7merm8 (N= ' + str(len(targets_7merm8)) +")"
    targets_7merm8.rename(columns={"7merm8": "type"}, inplace=True)
    targets_8mer['8mer'] = '8mer (N= ' + str(len(targets_8mer)) +")"
    targets_8mer.rename(columns={"8mer": "type"}, inplace=True)

    targets_6mer = targets_6mer[targets_6mer['b'].notna()]
    targets_7mera1 = targets_7mera1[targets_7mera1['b'].notna()]
    targets_7merm8 = targets_7merm8[targets_7merm8['b'].notna()]
    targets_8mer = targets_8mer[targets_8mer['b'].notna()]

    targets_6mer = targets_6mer.drop(targets_6mer.columns[[2]], axis=1)

    sleuth= sleuth[sleuth['b'].notna()]
    print(sleuth)
    sleuth = sleuth[~sleuth.target_id.isin(targets_6mer.target_id)]
    sleuth = sleuth[~sleuth.target_id.isin(targets_7mera1.target_id)]
    sleuth = sleuth[~sleuth.target_id.isin(targets_7merm8.target_id)]
    sleuth = sleuth[~sleuth.target_id.isin(targets_8mer.target_id)]
    print(sleuth)

    sleuth['se_b'] = 'non-targets (N= ' + str(len(sleuth)) +")"
    sleuth.rename(columns={"se_b": "type"}, inplace=True)

    final_df = pd.concat([targets_6mer, targets_7mera1, targets_7merm8, targets_8mer, sleuth])
    final_df = final_df.drop(final_df.columns[[0,3]], axis=1)
    print(final_df)

    final_df.to_csv(output, index=False)
