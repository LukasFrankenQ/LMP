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
    fix_zonal_remote_regions,
    update_daily,
    daily_func,
    prepare_system_total,
    prepare_household_total,
    prepare_constituency_total,
)
from _helpers import to_date_period, to_datetime
from _aggregation_helpers import flexible_aggregate


path = "results/periods/{}_{}.json"
template = "snakemake -call{} --configfile config/config.yaml -- {}"

target = "live/periods/{}_{}.json"

monthly_raw_path = "live/monthly_raw.json"
monthly_path = "live/monthly.json"

daily_raw_path = "live/daily_raw"
daily_path = "live/daily"
total_path = "live/total.json"

system_total_path = "live/system_total.json"
household_total_path = "live/household_total.json"
constituency_total_path = "live/constituency_total.json"

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

    now_time = to_datetime(day, period)

    daily_filename = str(Path(daily_raw_path) / f'{"-".join(day.split("-")[:-1])}.json')

    outfile = path.format(day, period)
    target = target.format(day, period)

    # os.system(template.format(" --touch", outfile))
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
                json.dump(daily, f, indent=4)

    except FileNotFoundError:
        shutil.copy(outfile, daily_filename)

    Path(daily_path).mkdir(parents=True, exist_ok=True)

    for fn in os.listdir(daily_raw_path):
        with open(Path(daily_raw_path) / fn, 'r') as f:        
            daily = json.load(f)

        daily = prepare_frontend_dict(daily, daily_func)
        fix_zonal_remote_regions(daily)

        with open(Path(daily_path) / fn, 'w') as f:
            json.dump(daily, f, indent=4)

    with open(monthly_raw_path, 'r') as f:
        monthly = json.load(f)

    monthly = update_monthly(new_step, monthly)
    fix_zonal_remote_regions(monthly)

    with open(monthly_raw_path, 'w') as f:
        json.dump(monthly, f, indent=4)

    print([pd.Timestamp.fromtimestamp(int(ts)) for ts in list(monthly)])
    regional_total = {list(monthly)[0]: flexible_aggregate(monthly)}

    system_total = prepare_system_total(
        regional_total,
        now_time,
        cfd_consumer_benefit_share=1.
        )

    with open(system_total_path, 'w') as f:
        json.dump(system_total, f, indent=4)

    household_total = prepare_household_total(
        regional_total,
        now_time,
        )
    
    with open(household_total_path, 'w') as f:
        json.dump(household_total, f, indent=4)

    const_total = prepare_constituency_total(
        household_total,
        'data/constituency_mapper.csv'
    )

    with open(constituency_total_path, 'w') as f:
        json.dump(const_total, f, indent=4)

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
