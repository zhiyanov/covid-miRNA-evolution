import sys
import pandas as pd


def cmp(df): 
    # print(df.columns)
    df["covid"] = df.index
    
    df["date"] = df["covid"].str.split("|").str[1]
    df["date"] = pd.to_datetime(df["date"], format="%Y-%m-%d")
    df["date"] = df["date"].astype("str").str.split(" ").str[0]
    df = df.sort_values(by=["date"])
    
    df = df.loc[df["date"].astype("str").apply(lambda dt: len(dt.split("-"))) == 3]
    
    first = df.drop(columns=["date", "covid"]).iloc[0].to_frame()
    first_date = df.iloc[0]["date"]
    first = first.rename(columns={first.columns.to_list()[0] : "positions"})
    first[first_date] = first["positions"].apply(lambda x: len(x.split(",")))
    
    last = df.drop(columns=["date", "covid"]).iloc[-1].to_frame()
    last_date = df.iloc[-1]["date"]
    last = last.rename(columns={last.columns.to_list()[0] : "positions"})
    last[last_date] = last["positions"].apply(lambda x: len(x.split(",")))

    result = first[[first_date]].join(last[[last_date]], how="inner")
    return result, first_date, last_date

def compare(*args):
    dfs = [pd.read_csv(path, sep=";", index_col=0) for path in args]
    print(args[0])
    result, first_date, last_date = cmp(dfs[0])

    for df, path in zip(dfs[1:], args[1:]):
        print(path)
        append, first_date, last_date = cmp(df) 
        protein = path.split("/")[-2].split("_")[-1]
        append = append.rename(columns={
            first_date: first_date + "_" + protein,
            last_date: last_date + "_" + protein,
        })
        result = result.join(append, how="inner")
    
    return result
