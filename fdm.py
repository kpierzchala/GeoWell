def d2tdr2(t, dx, a):
    out = (
        2.0
        * (
            (
                (t[1] - t[4]) * ((dx[0] * a[0] + dx[2] * a[1]))
                - (dx[1] * (t[4] - t[3]) * (dx[0] * a[3] + dx[2] * a[2]))
                / dx[3]
            )
        )
        / ((dx[0] + dx[2]) * ((dx[1]) ** 2.0 + dx[1] * dx[3]))
    )
    return out


def dtdr(t, dx, a, r):
    out = d2tdr2(t, dx, a) * r + (
        ((t[4] - t[3]) * (dx[0] * a[3] + dx[2] * a[2]))
        + ((dx[1] * dx[3]) / (dx[3] * dx[1] + (dx[1]) ** 2.0))
        * (
            (t[1] - t[4]) * (dx[0] * a[0] + dx[2] * a[1])
            + (t[3] - t[4]) * (dx[0] * a[3] + dx[2] * a[2])
        )
    ) / (dx[3] * (dx[0] + dx[2]))
    return out


def d2tdz2(t, dx, a, r):
    count1 = (
        (t[0] - t[4])
        * (
            ((r**2.0 - (r - dx[3]) ** 2.0) * a[3])
            + ((r + dx[1]) ** 2.0 - r**2.0) * a[0]
        )
        / ((r + dx[1]) ** 2.0 - (r - dx[3]) ** 2.0)
    )
    count2 = (
        (dx[0] / dx[2])
        * (t[2] - t[4])
        * (
            (
                (r**2.0 - (r - dx[3]) ** 2.0) * a[2]
                + ((r + dx[1]) ** 2.0 - r**2.0) * a[1]
            )
            / ((r + dx[1]) ** 2.0 - (r - dx[3]) ** 2.0)
        )
    )
    out = 2.0 * (count1 + count2) / ((dx[2]) ** 2.0 + dx[2] * dx[0])
    return out


def dtdtime(t, dx, a, r):
    out = d2tdr2(t, dx, a) + dtdr(t, dx, a, r) / r + d2tdz2(t, dx, a, r)
    return out
