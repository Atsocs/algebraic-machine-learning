def terms():
    alpha = ("B", "W")
    codes = [A + B + C + D for A in alpha for B in alpha for C in alpha for D in alpha]
    for code in codes:
        letters = tuple(map(lambda let: "white" if let == "W" else "black", code))
        A, B, C, D = letters
        text = "\\newcommand*{\\" + code + r"}{\term{\TermScale}"\
               + r"{1/1/" + B\
               + r", 1/2/" + A\
               + r", 2/1/" + C\
               + r", 2/2/" + D + r"}}"
        print(text)


if __name__ == '__main__':
    terms()
