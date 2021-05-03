#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Dieses \emph{Python}-Modul stellt alle Subroutinen bereit, die in direktem Zusammenhang mit dem betrachteten Problem stehen. Unter anderem verwendet es Matrixoperationen der \emph{scipy}-Bibliothek und das Datenbanksystem \emph{SQLite}.

import sqlite3
from math import sin, cos, pi, exp, sqrt
from scipy import zeros, empty, diag, linspace
from scipy.linalg import eigvalsh
from random import random


# \codesection{Modulweit gültige Konstanten}

# Die Systemparameter, die nach ihrer Initialisierung nicht mehr variiert werden, sind in globalen, d.h. hier modulweiten Variablen bzw. Konstanten gespeichert.

l = 10 # Größe der Superzelle

# \emph{on-site}-Energien (eV) für\dots
eC = 0.0 # Kohlenstoffatome
eX = 0.0 # Adsorbate

# \emph{hopping}-Parameter (eV) für\dots
t = 2.6 # benachbarte Kohlenstoffatome
V = 6.0 # Kohlenstoff-Adsorbat-Paare

# Nachbarschaftspenalty (eV)
p = 0.0


# \codesection{Orientierung im Wabengitter}

# Die Subroutinen \verb|ijk()| und \verb|xy()| vollziehen die Transformationen zwischen kartesischen Koordinaten \verb|(x, y)| und Index-Tripeln \verb|(i, j, k)| zur Verortung von Atomen. Die Modulo-Division in der Vorschrift für die $x$-Koordinate transformiert die rautenförmige in eine rechteckige Superzelle, deren Maße für $a = 1$ durch \verb|width()| und \verb|height()| gegeben sind.

def ijk(x, y, a = 1):
    x /= a
    y /= a
    
    i = (x / sqrt(3) - y / 3) % l
    j = y / 1.5 % l
    k = i % 1 + j % 1
    
    return int(i), int(j), int(k)

def xy(i, j, k, a = 1):
    x = a * sqrt(3) * ((i + 0.5 * (j + k)) % l + 0.5)
    y = a * (1.5 * j + 0.5 * k + 0.5)
    
    return x, y

width  = lambda: sqrt(3) * (l + 0.5)
height = lambda: 1.5 * l

# Atomkonfigurationen werden durch Listen \verb|X| von Tripeln \verb|(i, j, k)| dargestellt. \verb|neighbors()| weist jedem Atom seine nächsten drei Nachbarn zu. Sich so vortastend, gruppiert $\verb|cluster()|$ eine gegebene Konfiguration in Untermengen, die über nächste Nachbarschaft verbunden sind.

neighbors = lambda i, j, k: [
    ((i + 2 * k - 1) % l, j, 1 - k),
    (i, (j + 2 * k - 1) % l, 1 - k),
    (i, j, 1 - k)]

def cluster(X):
    free = X[:]
    tied = []
    while free:
        bunch = [free.pop()]
        for adsorbate in bunch:
            for neighbor in neighbors(*adsorbate):
                if neighbor in free:
                    bunch.append(neighbor)
                    free.remove(neighbor)
        tied.append(bunch)
    return tied

# \verb|carbons()| liefert die vollständige Konfiguration aller Kohlenstoffatome der Superzelle, sortiert von innen nach außen bezogen auf ihre Position in einem sechseckigen Cluster im darstellungsbezogenen Zentrum der Superzelle. So erhält man Cluster beliebiger Größe $\nX$, indem man einfach die ersten $\nX$ Elemente von \verb|carbons()| auswählt. Deswegen wird eine Atomkonfiguration hier auch als Liste und nicht als Menge (\verb|set|) von Tupeln aufgefasst, was logisch eigentlich zutreffender wäre, da die Reihenfolge physikalisch unbedeutend ist und keine doppelten Elemente auftreten dürfen. Das assoziative Array \verb|order| enthält ferner Sortierschlüssel um \verb|carbons()| so zu ordnen, dass zuerst die gewählte Standardkonfigurationen aufgefüllt wird.

def carbons():
    C = [(i, j, k)
        for i in range(l)
        for j in range(l)
        for k in 0, 1]

    i = C.index(ijk(0.5 * width(), 0.5 * height()))

    C[-1], C[i] = C[i], C[-1]
    
    return cluster(C)[0]

order = dict(
      C2X2 = lambda x: 1,
       C2X = lambda x: x[2],
       C6X = lambda x: (x[2], bool((x[1] - x[0]) % 3)),
    random = lambda x: random())


# \codesection{Berechnung der Energie}

# Der Großteil der eigentlichen Physik steckt im Tight-Binding-Modell, das schließlich auf den \textsc{Hamilton}-Operator bzw. die entsprechende Matrix führt. Letztere wird für eine gegebene Adsorbatkonfiguration von \verb|hamiltonian()| erzeugt. Ergänzend liefert \verb|energies()| direkt die nach Größe sortierten Energieeigenwerte, wobei jeder Wert doppelt vorkommt um der Entartung bezüglich des Spins gerecht zu werden. \verb|energy()| schließlich gibt für die Elektronenanzahl \verb|ne| nur die Gesamtenergie zurück, die bereits mit dem möglichen Penalty versehen ist. 

def hamiltonian(X = []):
    nC = 2 * l * l
    
    H = diag([eC] * nC + [eX] * len(X)) / 2
    
    for i in range(l):
        for j in range(l):
            a = 2 * (i * l + j)
            b = 2 * ((i + 1) % l * l + j)
            c = 2 * (i * l + (j + 1) % l)
            
            H[a + 1, (a, b, c)] = -t
    
    for n, (i, j, k) in enumerate(X):
        H[2 * (i * l + j) + k, nC + n] = V
    
    return H + H.T

def energies(X = []):
    return sorted(2 * list(eigvalsh(hamiltonian(X))))

def energy(X = [], ne = 0):
    return sum(energies(X)[:ne]) + penalty(X)

penalty = lambda X: p and p * 1.5 * K(X) * len(X)

# \verb|analytic()| berechnet die Eigenenergien des ungestörten Systems mit Hilfe der analytischen Lösung des Tight-Binding-Modells. Abgesehen von numerischen Ungenauigkeiten ist der Rückgabewert also identisch mit \verb|energies([])|.

def analytic():
    k = linspace(0, 2 * pi, l, endpoint = False)
    
    E = [t * sqrt(3 + 2 * (cos(k1) + cos(k2) + cos(k2 - k1)))
        for k1 in k for k2 in k]
    
    return sorted([eC + _ for _ in E] + [eC - _ for _ in E])


# \codesection{Bestimmung globaler Minima}

# \verb|combinations()| identifiziert iterativ alle möglichen Kombinationen von \verb|k| aus \verb|n| Elementen ohne Wiederholung in Form von $\tt n! \, / \, \bracks{k! \ (n - k)!}$ \verb|k|-Tupeln unterschiedlicher Indizes zwischen 0 und $\tt n - 1$ in aufsteigender Reihenfolge.

def combinations(N, K, n = 0, k = 1):
    if k > K:
        yield []
    else:
        for i in range(n, N - K + k):
            for _ in combinations(N, K, i + 1, k + 1):
                yield [i] + _

# \verb|best()| ermittelt, welche der vorgeschlagenen Konfigurationen die energetisch günstigste ist und gibt diese samt Spektrum und Gesamtenergie zurück. Diese Information \verb|I| kann bspw. auch an \verb|anneal()| weitergegeben werden.

def best(X, ne):
    I = { 'E' : float('inf') }

    for i, x in enumerate(X):
        J = { 'X' : x, 'e' : energies(x) }
        J['E'] = sum(J['e'][:ne]) + penalty(x)
        
        if J['E'] < I['E']:
            I = J
    
    return I


# \codesection{Simulated Annealing}

# Die Subroutine \verb|anneal()| startet ein Simulated Annealing, das nach \verb|n| in Folge abgelehnten bzw. maximal \verb|N| Schritten abgebrochen wird. \verb|R| ist die maximale Schrittweite, \verb|T0| die Ausgangstemperatur und \verb|f| die Kühlungsrate, die jeweils nach \verb|nf| Schritten Anwendung findet. Außerdem wird die Elektronenanzahl \verb|ne| erwartet. In \verb|I| sind die Konfiguration \verb|X|, die Energieeigenwerte \verb|e|, die Gesamtenergie \verb|E| sowie die Adsorbatbewegung vereint, wobei letztere durch eine Liste \verb|M| von Fremdatomtranspositionen und einen Index \verb|m| repräsentiert wird, der die aktuelle Position im Bewegungsablauf anzeigt. Die Werte werden durch \verb|anneal()| auf den Stand des optimalsten Zustands gebracht, der gefunden wurde, wobei erfolglose Folgeschritte der Vollständigkeit halber mit angehängt werden. Die Optimierung kann über \verb|on| auch ausgeschaltet werden. Dann findet nur die Initialisierung von \verb|I| inklusive der Berechnung der Energien statt, was manchmal nützlich sein kann. Ansonsten werden mit \verb|circle()| wiederholt alle Fremdatome durchlaufen.

def circle(n):
    while True:
        for i in xrange(n):
            yield i

def anneal(N = 1000, n = 10, R = 4.0, T0 = 0.0, f = 0.99, nf = 1, ne = 0,
    on = True, **I):
    
    if 'e' not in I: I['e'] = energies(I['X'])
    if 'E' not in I: I['E'] = sum(I['e'][:ne]) + penalty(I['X'])
    if 'M' not in I: I['M'] = []
    if 'm' not in I: I['m'] = len(I['M']) - 1
    
    if not (on and N and n and R > 0.75 and
        0 < len(I['X']) < 2 * l * l - 1): return I
    
    X  = I['X'][:]
    E0 = I['E']
    T  = T0
    
    i = 0
    j = 0
    
    for _ in circle(len(X)):
        r = random() * R
        a = random() * 2 * pi

        old = X[_]

        x, y = xy(*old)

        x += r * cos(a)
        y += r * sin(a)

        new = ijk(x, y)

        if new not in X:
            X[_] = new
            
            i += 1
            j += 1
            
            e = energies(X)
            E = sum(e[:ne]) + penalty(X)
        
            if E < E0 or T and random() < exp((E0 - E) / T):
                E0 = E
                j = 0
                
                I['M'].append([old, new])
            
                if E < I['E']:
                    I['e'] = e
                    I['E'] = E
                    I['X'] = X[:]
                    I['m'] = len(I['M']) - 1
            else:
                X[_] = old
            
            if i >= N or j >= n:
                break
            
            if not nf or not i % nf:
                T *= f
    return I


# \codesection{Größen zur Beschreibung der Adsorbatstruktur}

# Zur Berechnung der \textbf{Gruppiertheit} wird die Konfiguration zunächst mit \verb|cluster()| gruppiert.

def G(X = []):
    if len(X) < 2:
        return 0.0
    
    return             float(sum((len(x) - 1) * len(x)
        for x in cluster(X))) / ((len(X) - 1) * len(X))

# Um keine Kanten doppelt zu zählen wird die \textbf{Kompaktheit} nur von einem Untergitter ausgehend ermittelt.

def K(X):
    if not len(X):
        return 0.0
    
    K = 0.0
    
    for i, j, k in X:
        if k:
            for _ in neighbors(i, j, k):
                K += _ in X
    
    return 2.0 / 3.0 * K / len(X)

# Am aufwendigsten ist es, die \textbf{lokale Distanz} zu berechnen. Von jedem Fremdatom aus werden die jeweils umliegenden Kohlenstoffatome abgesucht, bis der nächste Nachbar gefunden ist. Solche gleicher "`Entfernung"' liegen dabei auf Sechsecken. Der Einfachheit halber wird hier vorrübergehend mit Mengen (\verb|set|) gerechnet.

def D(X):
    if not X:
        return 0.0
    
    if len(X) == 1:
        return 2.0 * l
    
    I = 0.0
    X = set(X)
    
    for x in X:
        O = [set(), {x}]
        
        while O[-1]:
            O.append({n for x in O[-1] for n in neighbors(*x)})
            
            O[-1] -= O[-3]
            
            if X & O[-1]:
                break
        
        I += len(O) - 2
    
    return I / len(X)

# Die Bestimmung der \textbf{Polarisierung} bzw. Untergitterasymmetrie ist hingegen sehr einfach, da nur die Kenntnis des Untergitterindex vonnöten ist.

def P(X = []):
    if not len(X):
        return 0.0
    
    return abs(1 - 2 * float(sum(k for i, j, k in X)) / len(X))

# \verb|describe()| liefert schließlich eine vollständige Beschreibung.

describe = lambda X: dict(G = G(X), K = K(X), D = D(X), P = P(X))


# \codesection{Berechnung der Zustandsdichte}

# Eine Darstellung der Zustanddichte wird durch \verb|dos()| generiert, indem jeder diskrete Eigenwert in \verb|E| durch eine \textsc{Lorentz}-Kurve der Halbwertsbreite \verb|dE| dargestellt wird. \verb|res| ist die Auflösung. Bei expliziter Angabe von \verb|ne| wird die Kurve an der \textsc{Fermi}-Grenze zweigeteilt.

def dos(E, dE = 0.05, res = 1001, ne = None):
    E = sorted(E)
    
    W = linspace(E[0], 2 * E[-1] - E[0], 2 * res - 1)
    delta =  1 / (dE ** 2 + (W - E[-1]) ** 2)
    f = (res - 1) / (E[-1] - E[0])
    
    delta *= f / (delta.sum() * len(E))
    
    DOS = zeros(res)
        
    for e in E:
        i = int(round(f * (E[-1] - e)))
        DOS += delta[i:res + i]
    
    I = dict(
          E = list(W[:res]),
        DOS = list(DOS))
    
    if ne is not None:
        mu = 0.5 * (E[max(ne, 1) - 1] + E[min(ne, len(E)) - 1])
        i = int(round(f * (mu - E[0])))
        
        I = [
            dict((k, v[:i + 1]) for k, v in I.items()),
            dict((k, v[    i:]) for k, v in I.items())]
    
    return I


# \codesection{Kovertierung}

# \verb|wX()| wandelt eine Konfiguration in Form einer Liste von Tripeln $(i \ j \ k)$ in eine Zeichenkette um, \verb|rX()| umgekehrt. \verb|wM()| und \verb|rM()| sind die Entsprechungen für Listen von Listen à zwei Tripeln \verb|(i, j, k)| und \verb|(I, J, k)|, wobei letztere je eine Fremdatomtransposition beschreiben. Diese Funktionen dienen der Kommunikation mit \emph{JavaScript}.

wX = lambda _: '~'.join('-'.join(str(_)
    for _ in _)
    for _ in _)

rX = lambda _: [tuple(int(_)
    for _ in _.split('-'))
    for _ in _.split('~')] if _ else []

wM = lambda _: ','.join('>'.join('-'.join(str(_)
    for _ in _)
    for _ in _)
    for _ in _)

rM = lambda _: [[tuple(int(_)
    for _ in _.split('-'))
    for _ in _.split('>')]
    for _ in _.split(',')] if _ else []


# \codesection{Scan des Phasendiagramms}

# Die Klasse \verb|Diagram| dient dem Zugriff auf die Bildpunkte des Phasendiagramms. Für einen gegebenen Bereich samt Auflösung stellt sie mit \verb|scan()| einen Iterator zur Verfügung.

class Diagram:
    def __init__(self,
        cX = (0.02, 0.4, 20),
        ce = (0.02, 0.4, 20)):
        
        self.nC = 2 * l * l
        
        self.size = (cX[2], ce[2])
        
        self.dce = max(ce[1] - ce[0], 1e-7) / (2 * ce[2] - 2 or 1)
        self.dcX = max(cX[1] - cX[0], 1e-7) / (2 * cX[2] - 2 or 1)
        
        self.cX = linspace(*cX)
        self.ce = linspace(*ce)
    
    def scan(self):
        for i, cX in enumerate(self.cX):
            nX = int(round(cX * self.nC))
        
            for j, ce in enumerate(self.ce):
                ne = int(round((1 + ce) * self.nC + nX))
            
                yield i, j, cX, ce, nX, ne


# \codesection{I/O}

# Den Abschluss bildet die Klasse \verb|IO| zur Kommunikation mit der \emph{SQLite}-Datenbank.

class IO:
    # Bei der Erzeugung einer Instanz \verb|io = IO()| über \verb|__init__()| wird die Verbindung zur Datenbank \verb|db| hergestellt, die dabei ggf. erst angelegt wird. Über \verb|io| kann dann nur auf eine Tabelle zugegriffen werden, die speziell für die zum Zeitpunkt der Initialisierung gültigen Systemparameter \verb|l|, \verb|eC|, \verb|t|, \verb|eX|, \verb|V| und \verb|p| angelegt wurde. Sie listet für verschiedene Kombinationen der Fremdatom- und Elektronenzahlen \verb|nX| und \verb|ne| die bislang günstigste Adsorbatkonfiguration \verb|X|, deren Gesamtenergie \verb|E| sowie die beschreibenden Größen \verb|G|, \verb|K|, \verb|D| und \verb|P| auf.
    
    def __init__(self, db = 'untitled.db'):
        self.db = db
        
        sqlite3.register_adapter(list, wX)
        sqlite3.register_converter('config', rX)
        
        self.connection = sqlite3.connect(db,
            detect_types = sqlite3.PARSE_DECLTYPES)
        
        self.connection.text_factory = str
        
        self.sql = self.connection.cursor()
        
        self.table = "'{l} {eC:g} {t:g} {eX:g} {V:g} {p:g}'".format(
                **globals())
        
        self.sql.execute('''
            create table if not exists {0} (
                nX integer, ne integer,
                 X config,   E float,
                 G float,    K float,
                 D float,    P float,
                primary key (nX, ne) on conflict replace
                ) without rowid'''.format(self.table))
    
    # \verb|get()| erwartet neben \verb|nX| und \verb|ne| weitere \emph{keyword arguments}. Die \emph{keys} geben die Namen der auszulesenden Spalten, die \emph{values} die Standardwerte an, für den Fall, dass die angeforte Zeile nicht existiert oder unvollständig ist.
    
    def get(self, nX, ne, **kwargs):
        if kwargs:
            keys = ', '.join(kwargs.keys())
        
            self.sql.execute(
                ''' select {0} from {1} where nX = ? and ne = ?
                '''.format(keys, self.table), [nX, ne])
        
            values = self.sql.fetchone()
        
            if values:
                kwargs = dict(zip(kwargs.keys(),
                    [b if b is not None else a
                    for a, b in zip(kwargs.values(), values)]))
        
        return kwargs
    
    # \verb|set()| legt nach einem analogen Schema Tabellenzeilen an bzw. ersetzt diese, falls sie bereits existieren. Die Werte der \verb|kwargs| werden in die jeweiligen Spalten geschrieben.
    
    def set(self, nX, ne, **kwargs):
        k = ', '.join(['nX', 'ne'] + kwargs.keys())
        v = ', '.join('?' * (2 + len(kwargs)))
        
        self.sql.execute(
            ''' insert into {0} ({1}) values ({2})
            '''.format(self.table, k, v), [nX, ne] + kwargs.values())
    
    # \verb|update()| ist wie \verb|set()|, nur dass hier ausschließlich die übergebenen Werte und nicht gleich die ganze Zeile überschrieben werden. Letztere muss dafür allerdings schon existieren.
    
    def update(self, nX, ne, **kwargs):
        do = ', '.join(_ + ' = ?' for _ in kwargs)
        
        self.sql.execute(
            ''' update {0} set {1} where nX = ? and ne = ?
            '''.format(self.table, do), kwargs.values() + [nX, ne])
    
    # \verb|map()| erwartet eine Instanz von \verb|Diagram|. Für jede der angeforderten Größen wird eine Matrix zurückgegeben, die die Werte für alle Bildpunkte des Phasendiagramms enthält.
    
    def map(self, diagram, **kwargs):
        A = dict()
        
        for key, value in kwargs.items():
            A[key] = empty(diagram.size, dtype = type(value))
            A[key].fill(value)
        
        for i, j, cX, ce, nX, ne in diagram.scan():                
            saved = self.get(nX = nX, ne = ne, **kwargs)
            
            for key, value in saved.items():
                if value is not None:
                    A[key][i, j] = value
        
        return A
    
    # Mit \verb|close()| werden schließlich Änderungen gespeichert und die Verbindung getrennt.
    
    def close(self):
        self.connection.commit()
        self.connection.close()