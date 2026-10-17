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
    if label == 'index.html':
        # 开局默认 = 128 体随机云（滑块也必须是 128，否则改了体数复现不出来）
        check(st['n'] == 128, 'index.html：开局默认 = 128 体', 'n=%s' % st['n'])
        gn = page.evaluate('() => ({v:+document.getElementById("genN").value,'
                         ' shown:document.getElementById("genNVal").textContent})')
        check(abs(gn['v'] - 2.10721) < 0.01 and gn['shown'] == '128',
              'index.html：#genN 滑块与默认场景一致', 'log10=%s 显示=%s' % (gn['v'], gn['shown']))
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
        # 预设 ↔ 表达式双向同步（手改表达式 → 预设跟随；不匹配 → 显示「自定义」）
        sy = pg.evaluate('''() => {
            const $ = id => document.getElementById(id);
            const setExpr = (id, v) => { const el=$(id); el.value=v; el.dispatchEvent(new Event('input')); };
            const pick = (id, v) => { const el=$(id); el.value=v; el.dispatchEvent(new Event('change')); };
            const o = {};
            o.initF = $('fPreset').value; o.initV = $('vPreset').value;
            setExpr('fExpr', 't^3');      o.matchF = $('fPreset').value;
            setExpr('vExpr', '1/r^3');    o.customV = $('vPreset').value;
            o.customShown = !$('vPreset').options[$('vPreset').options.length-1].hidden;
            setExpr('fExpr', '0.5*t+t^2'); o.customF = $('fPreset').value;
            pick('fPreset', 't');         o.backF = $('fExpr').value; o.backFdt = +$('dt').value;
            pick('vPreset', '1/r');       o.backV = $('vExpr').value;
            o.fAfter = $('fPreset').value; o.vAfter = $('vPreset').value;
            o.status = $('applyStatus').textContent;
            // 还原默认，后面的建云检查还在默认力律下做
            pick('fPreset', 't^2'); pick('vPreset', '1/r^2');
            return o;
        }''')
        check(sy['initF'] == 't^2' and sy['initV'] == '1/r^2',
              'index.html：启动时预设与表达式一致', 'f=%s F=%s' % (sy['initF'], sy['initV']))
        check(sy['matchF'] == 't^3', 'index.html：改表达式 → 预设跟随', 'f预设=%s' % sy['matchF'])
        check(sy['customV'] == '__custom__' and sy['customShown'],
              'index.html：非预设表达式 → 显示自定义', 'F预设=%s' % sy['customV'])
        check(sy['customF'] == '__custom__', 'index.html：f 同样支持自定义态', 'f预设=%s' % sy['customF'])
        check(sy['backF'] == 't' and sy['backV'] == '1/r',
              'index.html：选预设 → 回写表达式', 'f=%s F=%s' % (sy['backF'], sy['backV']))
        # dt 滑块 step=0.05，log10(2e-4)=-3.699 会被吸附到 -3.70，容差取半步
        check(abs(sy['backFdt'] - (-3.69897)) <= 0.025, 'index.html：选预设同时套用推荐 dt',
              'log10(dt)=%.4g' % sy['backFdt'])
        check('✓' in sy['status'], 'index.html：预设切换后内核重建成功', sy['status'])
        # 显示尺寸分级：N 越大越接近真实半径（真实半径必须与内核公式一致）
        sz = pg.evaluate('''() => {
            const d = +document.querySelector('#dimsSeg button.on').dataset.value;
            const pr = physRadiusFn();
            return {
                d: d,
                r1: pr(1), r8: pr(8),
                w: [sizeWeight(5), sizeWeight(1000), sizeWeight(2000),
                    sizeWeight(5000), sizeWeight(12000), sizeWeight(16384)],
                fl: [sizeFloor(100), sizeFloor(3000), sizeFloor(16384)],
                gl: [glowWeight(100), glowWeight(2000), glowWeight(5000)]
            };
        }''')
        import math
        rho = 100.0            # 密度滑块默认 10^2
        d = sz['d']
        exp1 = (1.0 * math.gamma(d / 2.0 + 1.0) / (rho * math.pi ** (d / 2.0))) ** (1.0 / d)
        exp8 = (8.0 * math.gamma(d / 2.0 + 1.0) / (rho * math.pi ** (d / 2.0))) ** (1.0 / d)
        check(abs(sz['r1'] - exp1) < 1e-9 and abs(sz['r8'] - exp8) < 1e-9,
              'index.html：真实半径 = (m·Γ(d/2+1)/(ρ·π^(d/2)))^(1/d)（与内核一致）',
              'd=%d R(1)=%.6g R(8)=%.6g' % (d, sz['r1'], sz['r8']))
        check(all(sz['w'][i] > sz['w'][i + 1] for i in range(len(sz['w']) - 1)) and sz['w'][0] == 1,
              'index.html：夸张权重随 N 单调下降（1 → 0.1）', str(sz['w']))
        check(sz['fl'][0] > sz['fl'][1] > sz['fl'][2] and sz['fl'][2] >= 1.0,
              'index.html：像素下限随 N 收缩但不小于 1px', str(sz['fl']))
        check(sz['gl'][0] > sz['gl'][1] > sz['gl'][2] and sz['gl'][2] == 0,
              'index.html：大 N 时光晕收掉', str(sz['gl']))
        pg.wait_for_timeout(400)
        pg.eval_on_selector('#genN', "el => { el.value = Math.log10(800); el.dispatchEvent(new Event('input')); }")
        pg.click('#genBtn')
        pg.wait_for_timeout(2500)
        st = pg.evaluate('() => ({n: Gravity.sim.n, G: Gravity.sim.G, e: Gravity.sim.e})')
        check(st['n'] > 100, 'index.html：随机云已生成', 'n=%s G=%.4g' % (st['n'], st['G']))
        check(st['e'] < 0, 'index.html：随机云是束缚态（E<0）', 'E=%s' % st['e'])
        # 滑块两端：下限 1、上限 16384（超 kMaxBodies 会静默丢体，只有数 n 才看得出）
        pg.eval_on_selector('#genN', "el => { el.value = 0; el.dispatchEvent(new Event('input')); }")
        pg.click('#genBtn')
        pg.wait_for_timeout(1500)
        n1 = pg.evaluate('() => Gravity.sim.n')
        check(n1 == 1, 'index.html：滑块下限 N=1 可建立', 'n=%s' % n1)
        pg.eval_on_selector('#genN', "el => { el.value = el.max; el.dispatchEvent(new Event('input')); }")
        pg.click('#genBtn')
        pg.wait_for_timeout(12000)
        st = pg.evaluate('() => ({n: Gravity.sim.n, len: Gravity.sim.bodies.length})')
        check(st['n'] == 16384, 'index.html：滑块上限 16384 体全部建成（含晕不丢体）', 'n=%s' % st['n'])
        check(st['len'] == st['n'] * 8, 'index.html：16384 体状态帧 = n×8', 'len=%s' % st['len'])
        # GPU 粒子层（WebGL2 点精灵）：要真的画出像素，且关掉后能回退到 Canvas 2D
        gl = pg.evaluate('''() => {
            const b = Gravity.sim.bodies, n = Math.floor(b.length / 8);
            const ok = gpu.ok;
            const drew = drawBodiesGPU(b, n, sizeWeight(n), sizeFloor(n), glowWeight(n), physRadiusFn());
            const glc = gpu.gl, W = glcv.width, H = glcv.height;
            const px = new Uint8Array(W * H * 4);
            glc.readPixels(0, 0, W, H, glc.RGBA, glc.UNSIGNED_BYTE, px);   // 必须在同一帧内读
            let cnt = 0; for (let i = 3; i < px.length; i += 4) if (px[i] > 8) cnt++;
            gpu.ok = false;                       // 回退路径必须返回 false（而不是抛异常）
            let fell = null;
            try { fell = drawBodiesGPU(b, n, sizeWeight(n), sizeFloor(n), glowWeight(n), physRadiusFn()); }
            catch (e) { fell = 'ERR ' + e.message; }
            gpu.ok = ok;
            return { ok: ok, drew: drew, painted: cnt, fallback: fell, maxPt: gpu.maxPt };
        }''')
        check(gl['ok'] is True and gl['drew'] is True, 'index.html：GPU 粒子层可用',
              'maxPt=%s' % gl['maxPt'])
        check(gl['painted'] > 1000, 'index.html：粒子层确实画出了像素',
              '非透明像素=%s' % gl['painted'])
        check(gl['fallback'] is False, 'index.html：无 WebGL 时回退 Canvas 2D（不抛异常）',
              str(gl['fallback']))
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
