if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Genetic screen workflow 0.0.1')
    parser.add_argument('--debug', dest='debug', required=False, type=int,
                        help='Debug level')
    parser.add_argument('--stage', dest='stage', required=False,
                        help='Limit job submission to a particular '
                             'Analysis stage. '
                             '[all,alignment,filter,duplicates,insertions,'
                             'annotate, grouping, count, plot]')
    parser.add_argument('--configuration', dest='configuration',
                        required=True, type=str, help='Configuration file (*.ini)')
    args = parser.parse_args()
