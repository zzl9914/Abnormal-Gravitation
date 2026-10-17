'use strict';
// 极简静态服务器（仿 Hash 2D-CAM）。
//
// 前端已经直连 wasm 内核，这里只负责把 web/ 托管出去——没有任何业务逻辑、
// 没有 SSE、没有 POST 接口。之所以还要一个服务器，是因为浏览器的
// WebAssembly / module 安全策略不允许 file:// 直接加载。
//
//   node web/tools/serve.js [--port 8123]

var http = require('http');
var fs = require('fs');
var path = require('path');
var url = require('url');

var webRoot = path.join(__dirname, '..');
var port = 8123;
var args = process.argv.slice(2);
for (var i = 0; i < args.length; i++) {
  if (args[i] === '--port') port = parseInt(args[++i], 10) || 8123;
}

var MIME = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.wasm': 'application/wasm',
  '.svg': 'image/svg+xml',
  '.ico': 'image/x-icon',
  '.png': 'image/png'
};

http.createServer(function (req, res) {
  var pathname = decodeURIComponent(url.parse(req.url).pathname || '/');
  if (pathname === '/' || pathname === '') pathname = '/index.html';
  // 防目录穿越
  var filePath = path.join(webRoot, path.normalize(pathname).replace(/^([/\\])+/, ''));
  if (filePath.indexOf(webRoot) !== 0) {
    res.writeHead(403);
    return res.end('forbidden');
  }
  fs.readFile(filePath, function (err, data) {
    if (err) {
      res.writeHead(404, { 'Content-Type': 'text/plain; charset=utf-8' });
      return res.end('404 ' + pathname);
    }
    res.writeHead(200, {
      'Content-Type': MIME[path.extname(filePath).toLowerCase()] || 'application/octet-stream',
      'Cache-Control': 'no-store',
      'Content-Length': data.length
    });
    res.end(data);
  });
}).listen(port, '127.0.0.1', function () {
  console.log('引力模拟已启动: http://127.0.0.1:' + port + '/   （Ctrl+C 停止）');
});
