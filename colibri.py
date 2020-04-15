
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument

if len(sys.argv) >= 2:
    config = sys.argv[1]
    vals = load_config(config)
    if vals:
        A, b, C, d, l, u = vals
        print("Okay")
        print(A, b, C, d, l, u)

else:
    m = 5
    A, b, C, d, l, u = generate(m)
    solver_lagrange_simple_b(m, A, b, C, d, l, u)
