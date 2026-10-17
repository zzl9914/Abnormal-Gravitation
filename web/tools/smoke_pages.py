# -*- coding: utf-8 -*-
"""两页 headless 冒烟：页面能加载 wasm 内核、能建模拟、能推进、无 JS 错误。

    python web/tools/smoke_pages.py            # 需要已启动 http://127.0.0.1:8123
"""
import sys
import time

from playwright.sync_api import sync_playwright

CHECKS = []


def check(cond, name, extra=''):
    CHECKS.append(bool(cond))
    print(('[PASS] ' if cond else '[FAIL] ') + name + ('  ' + extra if extra else ''))


def run(page, url, label, min_bodies):
    errors = []
    page.on('console', lambda m: errors.append(m.text) if m.type == 'error' else None)
    page.on('pageerror', lambda e: errors.append(str(e)))
    page.goto(url, wait_until='load')
    # 等待内核初始化 + 首次 setup
    page.wait_for_function('window.Gravity && window.Gravity.sim && window.Gravity.sim.n > 0',
                           timeout=15000)
    time.sleep(1.5)
    st = page.evaluate('''() => ({
        n: Gravity.sim.n, t: Gravity.sim.t, e: Gravity.sim.e, e0: Gravity.sim.e0,
        blen: Gravity.sim.bodies.length, running: Gravity.sim.running, G: Gravity.sim.G,
        err: document.getElementById('errbar') ? document.getElementById('errbar').textContent : ''
    })''')
    check(st['n'] >= min_bodies, '%s：天体已建立' % label, 'n=%s' % st['n'])
    check(st['blen'] == st['n'] * 8, '%s：状态帧 = n×8' % label, 'len=%s' % st['blen'])
    check(st['t'] > 0, '%s：时间在推进（内核在浏览器里跑）' % label, 't=%.4g' % st['t'])
    check(not st['err'], '%s：无错误提示条' % label, st['err'])
    # 渲染没有抛异常
    check(not errors, '%s：无 JS 错误' % label, '; '.join(errors[:3]))
    return st


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else 'http://127.0.0.1:8123'
    with sync_playwright() as p:
        # 优先用本机 Edge/Chrome（playwright 自带的 chromium 可能没下载）
        b = None
        for ch in ('msedge', 'chrome'):
            try:
                b = p.chromium.launch(channel=ch)
                break
            except Exception:
                b = None
        if b is None:
            b = p.chromium.launch()
        pg = b.new_page(viewport={'width': 1280, 'height': 800})
        run(pg, base + '/index.html', 'index.html', 2)
        # 主页：点「生成随机云」验证建云 + 标定链路（需要打开设置抽屉）
        pg.click('#settingsToggle')
        pg.wait_for_timeout(200)
        pg.eval_on_selector('#genN', "el => { el.value = Math.log10(800); el.dispatchEvent(new Event('input')); }")
        pg.click('#genBtn')
        pg.wait_for_timeout(2500)
        st = pg.evaluate('() => ({n: Gravity.sim.n, G: Gravity.sim.G, e: Gravity.sim.e})')
        check(st['n'] > 100, 'index.html：随机云已生成', 'n=%s G=%.4g' % (st['n'], st['G']))
        check(st['e'] < 0, 'index.html：随机云是束缚态（E<0）', 'E=%s' % st['e'])
        # 暂停 / 单步
        pg.click('#playBtn')
        pg.wait_for_timeout(300)
        check(pg.evaluate('() => Gravity.sim.running') is False, 'index.html：暂停生效')
        pg.click('#stepBtn')
        pg.wait_for_timeout(300)
        pg2 = b.new_page(viewport={'width': 1280, 'height': 800})
        run(pg2, base + '/three.html', 'three.html', 3)
        b.close()
    ok = all(CHECKS)
    print('\n' + ('ALL CHECKS PASSED' if ok else 'SOME CHECKS FAILED') +
          '  (%d passed, %d failed)' % (sum(CHECKS), len(CHECKS) - sum(CHECKS)))
    sys.exit(0 if ok else 1)


main()
