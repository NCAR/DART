import argparse
import datetime

def gregorian(date_pieces): 
    if len(date_pieces) == 2:
        gregorian_start = datetime.datetime(1601, 1, 1)
        delta = datetime.timedelta(
            days=int(date_pieces[1]),
            seconds=int(date_pieces[0])
        )
        ret_str = datetime.date.strftime(
            gregorian_start + delta,
            '%Y %-m %-d %-H %-M %-S'
        )
        ret_list = [int(ii) for ii in ret_str.split(' ')]
        return tuple(ret_list)

    if len(date_pieces) >= 3:
        delta = datetime.datetime(*date_pieces) - datetime.datetime(1601, 1, 1)
        return delta.seconds, delta.days


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Convert between %Y, %m,% d[, %H[, %M[, %S]]] and %S %d since 1601'
    )

    parser.add_argument(
        'date_pieces',
        metavar='',
        type=str,
        nargs='*'
    )
    args = parser.parse_args()
    date_pieces = [int(ii) for ii in args.date_pieces]    
    print(gregorian(date_pieces))
