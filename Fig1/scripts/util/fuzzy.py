import subprocess as process

class MatchResult:
    """scoring and storing helper for text matching results
    """
    def __init__(self, pos1: int, pos2: int,
                 val: float, text: str):
        self.pos1 = pos1
        self.pos2 = pos2
        self.val = val
        self.text = text
        self.match_length = len(text)

    def score(self):
        # length score(higher better) - proximity score(smaller better)
        return self.match_length - abs(self.pos1 - self.pos2)

    def iscloser(self, other: "MatchResult"):
        return abs(self.pos1 - self.pos2) < abs(other.pos1 - other.pos2)

    def islonger(self, other: "MatchResult"):
        return self.match_length > other.match_length

    def isbetter(self, other: "MatchResult"):
        return self.score() > other.score()


def fuzzy_phrase(text: str, regex: str, max_errors: int = 2) -> str | None:
    FUZZY_OPT = ["", f"Zbest{max_errors}"][max_errors > 0]
    try:
        out = process.check_output(f"echo \"{text}\" | ug -o{FUZZY_OPT} \"{regex}\"", shell=True)
        return out.decode("utf8").strip()
    except process.CalledProcessError as e:
        return None

def try_ugrep(text: str, regex: str) -> MatchResult | None:
    # number position, phrase position, number, phrase (separated by $)
    FORMAT = "--format=\"%[val]b$%[text]b$%[val]#$%[text]#%~\""
    try:
        out = process.check_output(f"echo \"{text}\" | ug -Po {FORMAT} \"{regex}\"", shell=True)
        out = out.decode("utf8").split("\n")[0].split(sep="$")
        return MatchResult(int(out[0]), int(out[1]), float(out[2]), str(out[3]))
    except process.CalledProcessError as e:
        return None


def take_best(text: str, current: MatchResult | None, pattern: str) -> MatchResult | None:
    result = try_ugrep(text, pattern)
    if (not current) or (current and result and result.isbetter(current)):
        return result
    else:
        return current


def classify_stage(text: str) -> tuple[float, str] | None:
    NUMBER_BEFORE = "(?<val>\\d+\\.?\\d*)[^0-9]*?"
    NUMBER_AFTER  = "[^0-9]*?(?<val>\\d+\\.?\\d*)"
    ret: MatchResult | None = None
    ret2: MatchResult | None = None
    phrase = fuzzy_phrase(text, "dpf|dfp", max_errors=0)
    if phrase:
        ret = take_best(text, ret, f"{NUMBER_BEFORE}(?<text>{phrase})")
        ret = take_best(text, ret, f"(?<text>{phrase}){NUMBER_AFTER}")

    phrase = fuzzy_phrase(text, "hpf|hfp", max_errors=0)
    if phrase:
        ret = take_best(text, ret, f"{NUMBER_BEFORE}(?<text>{phrase})")
        ret = take_best(text, ret, f"(?<text>{phrase}){NUMBER_AFTER}")

    phrase = fuzzy_phrase(text, "somites?", max_errors=2)
    if phrase:
        ret = take_best(text, ret, f"{NUMBER_BEFORE}(?<text>{phrase})")
        ret = take_best(text, ret, f"(?<text>{phrase}){NUMBER_AFTER}")

    phrase = fuzzy_phrase(text, "days? post[- ]?fertili[sz]ation", max_errors=2)
    if phrase:
        ret = take_best(text, ret, f"{NUMBER_BEFORE}(?<text>{phrase})")
        ret = take_best(text, ret, f"(?<text>{phrase}){NUMBER_AFTER}")

    phrase = fuzzy_phrase(text, "hours? post[- ]?fertili[sz]ation", max_errors=2)
    if phrase:
        ret = take_best(text, ret, f"{NUMBER_BEFORE}(?<text>{phrase})")
        ret = take_best(text, ret, f"(?<text>{phrase}){NUMBER_AFTER}")

    # note: less specific matches must come after more specific matches (hours post-fertilization)
    phrase = fuzzy_phrase(text, "hour", max_errors=1)
    if phrase:
        ret = take_best(text, ret, f"{NUMBER_BEFORE}(?<text>{phrase})")
        ret = take_best(text, ret, f"(?<text>{phrase}){NUMBER_AFTER}")

    if ret:
        return ret.val, ret.text
    else:
        return None