import os
import sys
import json
import shutil
import argparse
import pandas as pd
from pathlib import Path

sys.path.append(str(Path.cwd() / 'scripts'))
sys.path.append(str(Path.cwd() / '.github' / 'scripts'))

from _live_helpers import (
    update_monthly,
    prepare_frontend_dict,
    half_hourly_func,
    summary_func,
    fix_zonal_remote_regions,
    update_daily,
    daily_func,
)
from _helpers import to_date_period
from _aggregation_helpers import aggregate_stats


path = "results/periods/{}_{}.json"
template = "snakemake -call{} --configfile config/config.yaml -- {}"

target = "live/periods/{}_{}.json"

monthly_raw = "live/monthly_raw.json"
monthly_target = "live/monthly.json"

daily_raw_path = "live/daily_raw"
daily_target_path = "live/daily"
total_file = "live/total.json"
total_raw_file = "live/total_raw.json"

max_periods = 24


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="argparse to run day, periods that are not now"
        )

    parser.add_argument("--day", type=str, help="day to run", default=None)
    parser.add_argument("--period", type=int, help="period to run", default=None)

    args = parser.parse_args()
    day, period = args.day, args.period

    assert (int(day is None) + int(period is None)) % 2 == 0, "day and period must be both None or not None"

    if day is None:
        day, period = to_date_period(pd.Timestamp.now())

    daily_filename = str(Path(daily_raw_path) / f'{"-".join(day.split("-")[:-1])}.json')

    outfile = path.format(day, period)
    target = target.format(day, period)

    os.system(template.format(" --touch", outfile))
    os.system(template.format("", outfile))

    Path(target).parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(outfile, target)

    # load new data
    with open(target, 'r') as f:
        new_step = json.load(f)

    try:
        with open(daily_filename, 'r') as f:
            daily = json.load(f)
            daily = update_daily(daily, new_step, day)
            with open(daily_filename, 'w') as f:
                json.dump(daily, f)

    except FileNotFoundError:
        shutil.copy(outfile, daily_filename)

    Path(daily_target_path).mkdir(parents=True, exist_ok=True)

    for fn in os.listdir(daily_raw_path):
        with open(Path(daily_raw_path) / fn, 'r') as f:        
            daily = json.load(f)
        
        daily = prepare_frontend_dict(daily, daily_func)
        fix_zonal_remote_regions(daily)

        with open(Path(daily_target_path) / fn, 'w') as f:
            json.dump(daily, f)

    with open(monthly_raw, 'r') as f:
        monthly = json.load(f)

    monthly = update_monthly(new_step, monthly)
    # total = easy_aggregate(monthly)
    total = {list(monthly)[0]: aggregate_stats(monthly)}
    print('=====================================================')
    print('total raw')
    from pprint import pprint
    pprint(total[list(total)[0]]['national'])

    with open(monthly_raw, 'w') as f:
        json.dump(monthly, f)

    with open(total_raw_file, 'w') as f:
        json.dump(total, f)

    for data, func, fn in zip(
        [new_step, total, monthly],
        [half_hourly_func, summary_func, summary_func],
        [target, total_file, monthly_target]
        ):

        data = prepare_frontend_dict(data, func)
        fix_zonal_remote_regions(data)

        if 'total' in fn:
            # indicates last timestep at which total data was updated
            data[list(data)[0]]['last_update'] = list(new_step)[0]

        with open(fn, 'w') as f:
            json.dump(data, f)

    os.remove(outfile)

    # clean up periods, leaving only the last <max_periods> periods
    for fn in sorted(os.listdir("live/periods"))[:-max_periods]:
        os.remove(f"live/periods/{fn}")

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1, figsize=(10, 3.5))

    ax.text(0.5, 0.5, f"{day} {period}", fontsize=12, ha='center')

    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xticks([])
    ax.set_yticks([])

    plt.savefig('live/current_period.png')

    # add constituency totals
    with open(total_raw_file, 'r') as f:
        data = json.load(f)

    total_nodal = total[list(total)[0]]['nodal']['geographies']
    total_nodal = pd.DataFrame(
        {key: item['variables'] for key, item in total_nodal.items()}
        ).T.reset_index()
    total_nodal['index'] = total_nodal['index'].astype(str)

    total_zonal = total[list(total)[0]]['eso']['geographies']
    total_zonal = pd.DataFrame(
        {key: item['variables'] for key, item in total_zonal.items()}
        ).T
    
    const = pd.read_csv('data/constituency_mapper.csv')
    const['nodal_region'] = const['nodal_region'].astype(str)

    print(total_zonal.columns)
    const = const.merge(
        total_zonal[['single_rate_domestic']].reset_index(),
        left_on='zonal_region',
        right_on='index',
        how='left'
    ).drop(columns='index').rename(columns={'single_rate_domestic': 'single_rate_domestic_zonal'})

    const = (
            const.merge(
            total_nodal[['index', 'single_rate_domestic']],
            left_on='nodal_region',
            right_on='index',
            how='left'
        )
        .drop(columns='index')
        .rename(columns={'single_rate_domestic': 'single_rate_domestic_nodal'})
    )

    const_results = {}
    for l in ['nodal', 'zonal']:
        const[f'savings_{l}'] = const[f'single_rate_domestic_{l}'].mul(const['demand'])
        const_results[f'const_savings_{l}'] = const.groupby('constituency_code')[f'savings_{l}'].sum()

    f = pd.concat((
        const_results['const_savings_nodal'],
        const_results['const_savings_zonal'],
        ), axis=1)
        
    fdict = dict(f.T)

    for key, item in fdict.items():
        fdict[key] = {'variables': item.to_dict()}

    fdict = {list(total)[0]: {'constituencies': {'geographies': fdict}}}

    with open('live/constituency_total.json', 'w') as f:
        json.dump(fdict, f)
