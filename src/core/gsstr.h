#pragma once
// 最小字符/数字工具：freestanding 下没有 <cstring> / <cctype> / <cstdlib>，
// 自己实现词法与错误文本需要的一小撮函数。
#include <stddef.h>

inline bool gs_isspace(char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f';
}
inline bool gs_isdigit(char c) { return c >= '0' && c <= '9'; }
inline bool gs_isalpha(char c) {
    return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '_';
}
inline bool gs_isalnum(char c) { return gs_isalpha(c) || gs_isdigit(c); }

inline int gs_strlen(const char* s) {
    int n = 0;
    while (s[n]) n++;
    return n;
}

inline bool gs_streq(const char* a, const char* b) {
    while (*a && *a == *b) {
        a++;
        b++;
    }
    return *a == *b;
}

// 追加字符串到定长缓冲（永远以 '\0' 结尾，超长则截断）。返回当前长度。
inline int gs_strcat(char* dst, int cap, const char* src) {
    int n = 0;
    while (dst[n] && n < cap - 1) n++;
    for (int i = 0; src[i] && n < cap - 1; i++) dst[n++] = src[i];
    dst[n] = '\0';
    return n;
}

// 追加单个字符（定长缓冲，自动补 '\0'）
inline void gs_strcat_char(char* dst, int cap, char c) {
    int n = 0;
    while (dst[n] && n < cap - 1) n++;
    if (n < cap - 1) {
        dst[n] = c;
        dst[n + 1] = '\0';
    }
}

// 追加一个十进制数（最多 6 位有效数字，量级极小/极大时退化成定点），
// 用于错误信息里的 r / t 取值——不需要 printf 级别的精度。
inline int gs_strcat_num(char* dst, int cap, double x) {
    char buf[48];
    int n = 0;
    if (x < 0) {
        buf[n++] = '-';
        x = -x;
    }
    if (!(x == x)) return gs_strcat(dst, cap, "nan");
    long long ip = (long long)x;
    double frac = x - (double)ip;
    char tmp[24];
    int m = 0;
    if (ip == 0) tmp[m++] = '0';
    while (ip > 0 && m < 20) {
        tmp[m++] = (char)('0' + (int)(ip % 10));
        ip /= 10;
    }
    for (int i = m - 1; i >= 0 && n < 40; i--) buf[n++] = tmp[i];
    if (frac > 0.0 && n < 40) {
        buf[n++] = '.';
        for (int i = 0; i < 4 && n < 44; i++) {
            frac *= 10.0;
            int d = (int)frac;
            if (d > 9) d = 9;
            if (d < 0) d = 0;
            buf[n++] = (char)('0' + d);
            frac -= (double)d;
        }
    }
    buf[n] = '\0';
    return gs_strcat(dst, cap, buf);
}

// 极简 strtod：只支持十进制定点 + 可选指数（表达式字面量的全部需求）。
inline double gs_strtod(const char* s, const char** end) {
    const char* p = s;
    double sign = 1.0;
    if (*p == '+' || *p == '-') {
        if (*p == '-') sign = -1.0;
        p++;
    }
    double v = 0.0;
    bool any = false;
    while (gs_isdigit(*p)) {
        v = v * 10.0 + (double)(*p - '0');
        p++;
        any = true;
    }
    if (*p == '.') {
        p++;
        double scale = 0.1;
        while (gs_isdigit(*p)) {
            v += scale * (double)(*p - '0');
            scale *= 0.1;
            p++;
            any = true;
        }
    }
    if (!any) {
        if (end) *end = s;
        return 0.0;
    }
    if (*p == 'e' || *p == 'E') {
        const char* q = p + 1;
        bool eneg = false;
        if (*q == '+' || *q == '-') {
            eneg = (*q == '-');
            q++;
        }
        int e = 0;
        bool hasE = false;
        while (gs_isdigit(*q)) {
            e = e * 10 + (*q - '0');
            q++;
            hasE = true;
        }
        if (hasE) {
            p = q;
            double f = 1.0;
            for (int i = 0; i < e; i++) f *= 10.0;
            v = eneg ? v / f : v * f;
        }
    }
    if (end) *end = p;
    return sign * v;
}
