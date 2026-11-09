Backend API (`/predict`) behind

https://big-data-biology.org/software/macrel/

## Running

```bash
pixi run serve
```

This starts uwsgi on `127.0.0.1:8000` speaking the uwsgi protocol (not plain
HTTP), so it is meant to sit behind nginx (`uwsgi_pass 127.0.0.1:8000;`).
