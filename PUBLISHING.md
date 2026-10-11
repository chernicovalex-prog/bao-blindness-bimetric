# Публикация репозитория — инструкция (5 шагов)

Всё готово в папке `bao-blindness-bimetric/`. Осталась только твоя часть.

## 1. Создай репозиторий на GitHub
- github.com → New repository → имя `bao-blindness-bimetric`
- Public, БЕЗ галочек (не добавляй README/license — они уже есть)

## 2. Запушь код
```bash
cd bao-blindness-bimetric
git init
git add -A
git commit -m "v10.6: canonical pipeline (BAO/SNe/CMB) + CI + verification artefacts"
git branch -M main
git remote add origin git@github.com:<ТВОЙ_ЛОГИН>/bao-blindness-bimetric.git
git push -u origin main
```
(нужен SSH-ключ в настройках GitHub, либо https + token)

## 3. Проверь CI
- Вкладка Actions репозитория: через пару минут оба job'а (checks, full-scan)
  должны стать зелёными. Если что-то упадёт — пришли мне лог.

## 4. Zenodo (DOI)
- zenodo.org → Log in with GitHub → репозиторий в списке → переключатель ON
- На GitHub: Releases → Draft a new release → тег `v10.6` → Publish
- Zenodo автоматически создаст DOI

## 5. Пришли мне DOI
- Я вставлю его в статью (v21.6, Data Availability + CITATION) одной
  транзакцией с Журналом/Дневником.

## Примечания
- Данные Pantheon+ в репо не включены (лицензия); `make data` скачивает их
  с официального GitHub-репозитория Pantheon+SH0ES. CI их не тянет — SNe-прогон
  остаётся локальным (`make sne-full`).
- Слепые верификационные коды лежат в `verification/` с хэшем и логами.
