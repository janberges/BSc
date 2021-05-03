#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Dieses \emph{CGI}-Skript regelt die Kommunikation mit dem Benutzer. Es greift dafür auf Subroutinen aus \emph{graphene.py} und \emph{graphics.py} zurück.

import cgi;
form = cgi.FieldStorage()

import cgitb;
cgitb.enable()

import graphene as gn
import graphics as gx

import os
import re

names = globals()


# \codesection{Globale Parameter}

# Alle Parameter werden jeweils dort eingelesen, wo sie verarbeitet werden. Die Schnittstelle liefert die Inhalte der Formularelemente des HTML-Dokuments als Zeichenketten, die noch in die jeweiligen Typen konvertiert werden müssen. Dabei werden Ganzzahlen mit \verb|int(float())| statt mit \verb|int()| verarbeitet, um einen Fehler bei der Eingabe einer Kommazahl wie "`1.0"' zu verhindern. Zunächst werden hier die globalen Systemkonstanten (siehe \emph{graphene.py}) angefordert, wobei teilweise Einschränkungen des Wertebereichs vorgenommen werden.

l = gn.l = min(max(int(float(form.getfirst('l', 12))), 1), 30)

nC = 2 * l * l

eC = gn.eC = float(form.getfirst('eC', 0.0))
eX = gn.eX = float(form.getfirst('eX', 0.0))
t  = gn.t  = float(form.getfirst( 't', 2.6))
V  = gn.V  = float(form.getfirst( 'V', 6.0))
p  = gn.p  = float(form.getfirst( 'p', 0.0))

nX = min(max(int(float(form.getfirst('nX',       1))), 0),            nC)
ne = min(max(int(float(form.getfirst('ne', nC + nX))), 0), 2 * (nC + nX))


# \codesection{Verbindung zur Datenbank}

# Nachdem eine Verbindung zur \emph{SQLite}-Datenbank aufgebaut wurde, dürfen die globalen Variablen nicht mehr verändert werden. Ansonsten würden die Informationen zu einem gefundenen Minimum in eine falsche Tabelle geschrieben werden.

db = str(form.getfirst('db', 'example.db'))

if not os.path.exists(db):
    db = 'empty.db'

io = gn.IO(db)

# Um stets zu wissen, ob eine bessere Konfiguration gefunden wurde, werden die Daten des bisherigen Optimums geladen. Sind diese noch nicht vorhanden, wird das erreichte Energieoptimum auf $\infty$ gesetzt.

old = io.get(nX = nX, ne = ne, E = float('inf'), X = [])


# \codesection{Die aktuelle Adsorbatkonfiguration}

# Jetzt wird die vom Benutzer übermittelte Adsorbatkonfiguration \verb|X| eingelesen und angepasst, falls sie nicht vorhanden, unvollständig, fehlerhaft oder überzählig ist.

X = gn.rX(form.getfirst('X', ''))

# Die Positionen aller Kohlenstoffatome werden in \verb|C| gespeichert und nach einem Vergleich mit \verb|X| doppelte Fremdatome oder solche, die sich außerhalb der Superzelle befinden entfernt.

C = gn.carbons()
X = list(set(_ for _ in X if _ in C))

# Falls keine Konfiguration übermittelt wurde, kommt die gewählte Voreinstellung zum tragen: Soll, wenn vorhanden, die zuletzt abgespeicherte Konfiguration angezeigt werden (\verb|saved|) oder die energetisch günstigste der ausgewählten Standardkonfigurationen \verb|X0|?

saved = ' checked' if 'saved' in form or not form else ''
X0 = form.getlist('X0') or ['C2X2']

if not X:
    if old['X'] and saved:
        names.update(old)
    else:
        selection = [sorted(C, key = gn.order[_])[:nX] for _ in X0]
        names.update(gn.best(selection, ne = ne))

# Sind für gegebenes \verb|nX| hingegen zu viele oder zu wenige Fremdatome übermittelt worden, wird die entsprechende Anzahl entfernt bzw. mittig hinzugefügt.

elif len(X) > nX: X[nX:] = []
elif len(X) < nX: X.extend([_ for _ in C if _ not in X][:nX - len(X)])

# Nur wenn eine eingegebene Konfiguration unverändert übernommen wurde, macht es Sinn auch den Bewegungsverlauf in Form von \verb|M| und \verb|m| zu laden.

else:
    M =     gn.rM(form.getfirst('M', ''))
    m = int(float(form.getfirst('m', -1)))


# \codesection{Simulated Annealing}

N  = min(max(int(float(form.getfirst( 'N', 1000))),  0), 1000)
n  =     max(int(float(form.getfirst( 'n',  100))),  0)
R  =     max(    float(form.getfirst( 'R',  4.0)), 1.0)
T0 =             float(form.getfirst('T0',  0.0))
f  =             float(form.getfirst( 'f', 0.99))
nf =     max(int(float(form.getfirst('nf',    1))),  1)

# Im Anschluss ans Simulated Annealing werden \verb|X|, \verb|e|, \verb|E|, \verb|M| und \verb|m| an den besten Zustand angepasst, der gefunden wurde. Mögliche erfolglose Folgeschritte sind in \verb|M| mit aufgelistet, befinden allerdings von \verb|m| aus gesehen in der "`Zukunft"' und können ignoriert werden.

on = True if 'anneal' in form else False

names.update(gn.anneal(**names))


# \codesection{Beschreibung und Einpflegen möglicher Fortschritte}

# Die neue Konfiguration \verb|X| kann nun beschrieben und, falls sie die bislang energetisch günstigste ist, in die Datenbank eingelesen werden. Ein Speichern erfolgt hier noch nicht.

how = gn.describe(X)

if E < old['E'] - 1e-7:
    io.set(nX = nX, ne = ne, X = X, E = E, **how)


# \codesection{Darstellung der Superzelle}

# Hier wird ein Teil der \emph{SVG}-Grafik erstellt, die später die Superzelle darstellt.

width = 450
a = width / gn.width()
height = a * gn.height()

space = '\n' + ' ' * 16

c = "<circle fill='url(#grey)' cx='{0}' cy='{1}' r='%g' />" % (a / 4)
x = "<circle fill='url(#blue)' cx='{0}' cy='{1}' r='%g'" % (.3 * a) \
    + " id='{2}-{3}-{4}' onmousedown='grab(this); return false' />"

C = space.join(c.format(*gn.xy(*_, a = a)    ) for _ in C)
X = space.join(x.format(*gn.xy(*_, a = a) + _) for _ in X)


# \codesection{Allgemeine Ausgabeparameter}

# Die Größe und das Aussehen der Abbildungen werden an die restliche Benutzeroberfläche angepasst. Hier ist zudem eine weitere Formatierungsvorschrift für das Abspeichern von Abbildungen angegeben. Letzteres wird durch \verb|savefigures| ein- und ausgeschaltet.

screen = dict(
    font_family = 'Verdana, Geneva, sans-serif',
    font_size = 12.0,
    stroke_width = 0.7,
    margin_left = 50)

textwidth = 398.33861 # pt
halfwidth = 0.48 * textwidth

latex = dict(
    save = True,
    font_family = 'Iwona',
    font_size = 9.0,
    width = halfwidth,
    height = 2.0 / 3 * halfwidth,
    stroke_width = 0.5,
    x_tick_spacing = 40,
    y_tick_spacing = 40,
    margin_left = None)

savefigures = False


# \codesection{Darstellung der Zustandsdichte}

# Die Eigenenergien \verb|e| der aktuellen Adsorbatkonfiguration werden in Einheiten von \verb|t| umgerechnet und als Zustandsdichte dargestellt. Optional ist zum Vergleich die simultane Darstellung des Spektrums reines Graphens möglich.

e = [_ / t for _ in e]

dE  =     max(    float(form.getfirst( 'dE', 0.02)), 1e-7)
res = min(max(int(float(form.getfirst('res', 1001))),   3), 3000)

dos = gn.dos(E = e, dE = dE, res = res, ne = ne)

# Soll auch die DOS reinen Graphens dargestellt werden?

clean = ' checked' if 'clean' in form or not form else ''

if clean:
    dos.append(
        gn.dos(E = [_ / t for _ in gn.analytic()], dE = dE, res = res))

legend = ['occupied', 'vacant', 'clean graphene']
color  = [    '#00f',   '#888',           '#f80']

plots = [
    dict(
        x = dos[_]['E'],
        y = dos[_]['DOS'],
        legend = legend[_],
        line = dict(
            stroke = color[_],
            fill = color[_],
            fill_opacity = 0.5))
    for _ in reversed(range(len(dos)))]

figure = dict(
    width = 550,
    height = 300,
    x_label = '[E] / [t]',
    y_label = '[DOS]',
    y_ref = 0,
    plots = plots,
    **screen)

dos = gx.plot(**figure)

# Im Folgenden wird die optionale Erzeugung einer SVG-Datei eingeleitet:

if savefigures:
    figure.update(latex)
    gx.plot(name = 'fig/DOS', **figure)


# \codesection{Darstellung des Phasendiagramms}

# Zunächst wird eingelesen, welche Größen durch welche Farbe für welchen Bereich des Phasendiagramms bei welcher Auflösung dargestellt werden sollen.

channels = dict(
    red   = str(form.getfirst(  'red', 'K')),
    green = str(form.getfirst('green', 'G')),
    blue  = str(form.getfirst( 'blue', 'D')))

ce0 =             float(form.getfirst('ce0', 0.02))
cX0 =             float(form.getfirst('cX0', 0.02))
ce1 =             float(form.getfirst('ce1',  0.4))
cX1 =             float(form.getfirst('cX1',  0.4))
nce = min(max(int(float(form.getfirst('nce',   20))), 1), 100)
ncX = min(max(int(float(form.getfirst('ncX',   20))), 1), 100)

# \verb|scan| gibt an, ob momentan ein Scan des Phasendiagamms stattfindet. \verb|goto| identifiziert für diesen Fall den Pixel, der als nächster prozessiert werden soll.

scan = str(form.getfirst('scan', ''))
goto = ''

# In den folgenden Variablen werden die Abbildungen des Phasendiagramms und möglicher Querschnitte durch dieses gespeichert.

phases = ''
cut = { 'row' : '', 'col' : '' }

# Sollen diese Abbildungen überhaupt erzeugt werden?

showphases = ' checked' if 'showphases' in form or not form else ''

if showphases:
    # Zunächst wird eine Instanz von \verb|Diagram| erzeugt, die neben dem Iterator auch die Werte bereitstellt, die auf die Koordinatenachsen aufgetragen werden.

    diagram = gn.Diagram(
        cX = (cX0, cX1, ncX),
        ce = (ce0, ce1, nce))

    # \verb|show| speichert die Namen und Standardwerte der darzustellenden Größen, deren Werte in Matrizen eingelesen werden, die in \verb|A| vereint sind.

    show = dict((v, 0.0) for v in channels.values() if v != 'none')

    A = io.map(diagram, **show)

    # Zur Farbgebung müssen Minima und \emph{peak-to-peak}-Werte aller Größen bekannt sein.

    inf = dict((_, A[_].min()              ) for _ in show)
    ptp = dict((_, A[_].max() - inf[_] or 1) for _ in show)

    # In \verb|pixels| werden die Daten zur Darstellung jedes Pixel gesammelt. Wenn es einen Pixel gibt, der den aktuellen Werten von \verb|nX| und \verb|nE| entspricht, speichert \verb|selected| dessen Indizes. Sein Auffinden veranlasst zudem das Warten auf den folgenden Pixel, der als nächster gescannt würde.

    pixels   = []
    selected = None
    wait     = False

    # Hier beginnt die Schleife, in der das Phasendiagramm abgelaufen wird.

    for i, j, cX, ce, mX, me in diagram.scan():
        # Als erstes werden die Intensitäten für den roten, grünen und blauen Farbkanal berechnet. 0 entspricht dem Minimum, 255 dem Maximum der jeweiligen Größe, unabhängig vom \emph{peak-to-peak}-Wert.

        rgb = dict(
            (k, int(255 * (A[v][i, j] - inf[v]) / ptp[v])
            if v in A else 0) for k, v in channels.items())

        rgb = 'rgb({red}, {green}, {blue})'.format(**rgb)

        # Jeder Pixel erhält eine \verb|ID|, die sich aus \verb|nX| und \verb|nE| zusammensetzt. Ist er nicht \verb|selected|, gehört er der \emph{CSS}-Klasse \verb|px| an.

        ID = '{}-{}'.format(me, mX)
        Class = 'px'

        # Nach dem Warten wird der nächste Pixel identifiziert\dots

        if wait:
            goto = ID
            wait = False

        # \dots{}wobei ersteres durch das Auffinden des ausgewählten Pixels getriggert wurde.

        if (ne, nX) == (me, mX):
            selected = i, j
            Class = 'selected'
            if scan:
                wait = True

        # Im \emph{Tooltip} sollen die Werte aller ausgewählten Größen stehen.

        title = '&#10;'.join('{} = {:.4g}'.format(v, A[v][i, j])
            for k, v in channels.items() if v in A)

        pixels.append(
            dict(
                 x = ce,
                 y = cX,
                dx = diagram.dce,
                dy = diagram.dcX,
                bar = dict(
                    Class = Class,
                    id = ID,
                    fill = rgb,
                    stroke = rgb,
                    onclick = 'lookat(this.id)',
                    title = title)))

    figure = dict(
        width = 550,
        height = 400,
        x_label = '[c]_{e}',
        y_label = '[c]_{X}',
        plots = pixels,
        **screen)

    phases = gx.plot(**figure)

    if savefigures:
        figure.update(latex)
        gx.plot(name = 'fig/phases', **figure)


    # \codesection{Darstellung der Querschnitte}

    # Wenn ein Pixel gewählt ist, werden "`Höhenprofile"' der Zeile und der Spalte des Phasendiagramms, denen er angehört, erzeugt.

    if selected and show:
        i, j = selected

        # Auf der Ordinate sind dabei alle ausgewählten Größen dargestellt.

        ylabel = ' &#8239; '.join('<{0}>{{{1}}}'.format(k, v)
            for k, v in channels.items() if v != 'none')

        for name, xlabel, x, indices in [
            ('row', '[c]_{e}', diagram.ce, (i, slice(None))),
            ('col', '[c]_{X}', diagram.cX, (slice(None), j))]:

            figure = dict(
                width = 275,
                height = 180,
                x_tick_spacing = 55,
                y_tick_spacing = 45,
                x_label = xlabel,
                y_label = ylabel,
                y_ref = 0,
                plots = [
                    dict(
                        x = x,
                        y = A[v][indices],
                        line = dict(
                            fill_opacity = 0.2,
                            stroke = k,
                            fill = k))
                    for k, v in channels.items() if v in A],
                **screen)

            cut[name] = gx.plot(**figure)

            if savefigures:
                figure.update(latex)
                figure['y_tick_spacing'] = 30
                gx.plot(name = 'fig/' + name, **figure)


# \codesection{Speichern und Trennung von der Datenbank}

# Die Verbindung zur Datenbank kann jetzt getrennt werden. Mögliche Fortschritte werden nur dann gespeichert, wenn der Schreibzugriff auf die Datenbank explizit erlaubt wurde.

write = ' checked' if 'write' in form else ''

if write:
    io.close()


# \codesection{Formatierung der Ausgabe}

# Schließlich werden noch die Gesamtenergie \verb|E| und der Bewegungsverlauf \verb|M| formatiert sowie die Auswahlmenüs für die Farbkanäle und die Checkbox-Selektoren generiert.

E = '{:g}'.format(E).replace('-', '&#8722;')
M = gn.wM(M)

for k, v in channels.items():
    channels[k] = '''
      <select name='{}'>
          <option value='none'>-</option>
          <option value='E'>total energy</option>
          <option value='D'>local distance</option>
          <option value='K'>compactness</option>
          <option value='G'>aggregation</option>
          <option value='P'>polarisation</option>
      </select>
     '''.format(k).replace(
         "'{}'".format(v), "'{}' selected".format(v))

X0 = dict((key, ' checked' if key in X0 else '')
    for key in gn.order.keys())

# Jetzt können alle Lücken im HTML-Template gefüllt und dieses an den Browser gesendet werden. Für den Fall, dass auch eine SVG-Datei der Superzelle erzeugt werden soll, muss ein kleiner Umweg gemacht werden.

print 'Content-type: text/html\n'

with open('structure.html') as template:
    if savefigures:
        html = template.read().format(**vars())

        print html

        with open('fig/supercell.svg', 'w') as svg:
            svg.write(re.search('(<svg.+?/svg>)', html, re.S).group())
    else:
        print template.read().format(**vars())
