import os
import sys
import json
import shutil
import pandas as pd
from pathlib import Path

sys.path.append(str(Path.cwd() / 'scripts'))
sys.path.append(str(Path.cwd() / '.github' / 'scripts'))

from _live_helpers import (
    update_monthly,
    prepare_frontend_dict,
    half_hourly_func,
    summary_func,
    easy_aggregate,
)
from _helpers import to_date_period


path = "results/periods/{}_{}.json"
template = "snakemake -call{} --configfile config/config.yaml -- {}"

target = "live/periods/{}_{}.json"
monthly_file = "live/monthly.json"
total_file = "live/total.json"

max_periods = 24


if __name__ == "__main__":

    now = pd.Timestamp.now()
    day, period = to_date_period(now)

    outfile = path.format(day, period)
    target = target.format(day, period)

    os.system(template.format(" --touch", outfile))
    os.system(template.format("", outfile))

    Path(target).parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(outfile, target)

    # load new data
    with open(target, 'r') as f:
        new_step = json.load(f)

    with open(monthly_file, 'r') as f:
        monthly = json.load(f)

    monthly = update_monthly(new_step, monthly)
    total = easy_aggregate(monthly)

    for data, func, fn in zip(
        [new_step, total, monthly],
        [half_hourly_func, summary_func, summary_func],
        [target, total_file, monthly_file]
        ):

        data = prepare_frontend_dict(data, func)

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
    plt.show()
