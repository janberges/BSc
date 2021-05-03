// Das foglende \emph{JavaScript} ist für die dynamischen Elemente und Animationen auf der Benutzeroberfläche verantwortlich. Zudem greift es lesend und schreibend auf das HTML-Formular zu, das der Kommunikation mit \emph{index.py} dient.


// \codesection{Globale Variablen und \emph{Event-Handler}}

atom = null  // manuell bewegtes Atom
scan = false // Scan läuft
play = false // Film läuft

step = new Array() // aktive Animationen nach Atom-ID

X = Y = 0     // Koordinaten des Cursors im Browserfenster
I = J = K = 0 // manuelles Bewegen: Indizes der \textsc{Voronoi}-Region
dx = dy = 0   // Differenz zwischen Cursor und Atomkoordinaten

document.onmousemove = hold // Atom festhalten
document.onkeypress  = type // Tastatureingabe

// Nach dem Laden des HTML-Inhalts werden ausgewählte Elemente ausgelesen.

onload = function() {
    // Zunächst werden zwei Abkürzungen definiert.

    form = document.form                   // Formular
    info = document.getElementById('info') // $(i \ j \ k)$-Anzeige

    l =   parseInt(form.l.value) // Größe der Superzelle
    a = parseFloat(form.a.value) // Atomabstand
    M =         rM(form.M.value) // Bewegungsverlauf/Transpositionen
    m =   parseInt(form.m.value) // aktueller Index von \verb|M|

    // Gegebenenfalls wird der Scan des nächsten Bildpunkts eingeleitet.

    if (form.scan.value && form.goto.value) {
        form.anneal.value = 1
        scan = setTimeout(function() { lookat(form.goto.value) }, 500)
        }
    }


// \codesection{Methoden für Zahlen}

// Der implementierte Modulo-Operator \verb|%| von JavaScript ist für negative Zahlen fehlerhaft. Diese Methode (von \textsc{Stephen Chapman}, \url{http://javascript.about.com/od/problemsolving/a/modulobug.htm}, 28. August 2014) schafft Abhilfe:

Number.prototype.mod = function(n) {
    return ((this % n) + n) % n
    }

// \verb|x.toFixed(n)| wandelt die Zahl \verb|x| in eine Zeichenkette mit n Nachkommastellen um. Um mögliche angehängte Nullen zu entfernen, wird das Ergebnis wieder als Zahl interpretiert.

Number.prototype.to = function(n) {
    return parseFloat(this.toFixed(n))
    }


// \codesection{Konvertierung}

// \verb|wM()| wandelt die Liste der Transpositionen in eine Zeichenkette um, \verb|rM()| umgekehrt.

function wM(x) {
    return x.map(function(y) { return y.join('>') }).join(',')
    }

function rM(x) {
    return x ? x.split(',').map(function(y) { return y.split('>') }) : []
    }

// \verb|getX()| und \verb|getY()| liefern kartesische Koordinaten für ein gegebenes Index-Tripel \verb|[i, j, k]|.

function getX(i, j, k) {
    return a * Math.sqrt(3) * ((i + .5 * (j + k)).mod(l) + .5)
    }

function getY(i, j, k) {
    return a * (1.5 * j + .5 * k + .5)
    }


// \codesection{Absenden des Formulars}

// Vor dem Senden müssen zunächst die Daten zur Adsorbatkonfiguration und deren Verlauf in Textform gebracht werden, da diese sowohl von JavaScript als auch Python verstanden wird.

function send() {
    var atoms = document.getElementsByTagName('circle');
    var ids = new Array()

    for (var i = 0; i < atoms.length; i++)
        if (atoms[i].id)
            ids.push(atoms[i].id)

    form.X.value = ids.join('~')
    form.M.value = wM(M)
    form.m.value = m
    form.submit()
    }

// \verb|lookat(id)| lädt die Daten für einen Bildpunkt gegebener \verb|id|. Informationen zur aktuellen Adsorbatstruktur sowie Änderungen des Formularinhalts werden dabei verworfen.

function lookat(id) {
    var n = id.split('-')

    form.l.value  = form.l.defaultValue
    form.eC.value = form.eC.defaultValue
    form.eX.value = form.eX.defaultValue
    form.t.value  = form.t.defaultValue
    form.V.value  = form.V.defaultValue

    form.ne.value = n[0] // Die \verb|id| selbst liefert \verb|nE|\dots
    form.nX.value = n[1] // \dots{}und \verb|\nX|

    form.submit()
    }


// \codesection{Manuelle Fremdatombewegung}

// Mit \verb|grab()| wird ein Objekt, hier ein Fremdatom, "`aufgehoben"'. Die Position des übergeordneten \emph{SVG}-Elements ist egal, da die Differenz zwischen Cursor- und Objektkoordinaten gespeichert wird.

function grab(object) {
    atom = object

    dx = X - atom.getAttribute('cx')
    dy = Y - atom.getAttribute('cy')

    document.onmouseup = drop
    }

// \verb|hold()| bestimmt permanent die aktuelle Cursorposition. Wird zudem gerade ein Fremdatom "`festgehalten"', wird dessen Position angepasst. Gleichzeitig werden dann die Indizes der \textsc{Voronoi}-Region, über der sich das Atom gerade befindet, neben dem Cursor angezeigt.

function hold(event) {
    event = event || window.event

    X = event.clientX
    Y = event.clientY

    if (atom) {
        x = X - dx
        y = Y - dy

        i = (x / Math.sqrt(3) - y / 3) / a
        j = y / 1.5 / a

        I = Math.floor(i.mod(l))
        J = Math.floor(j.mod(l))
        K = Math.floor(i.mod(1) + j.mod(1))

        atom.setAttribute('cx', x)
        atom.setAttribute('cy', y)

        marker.setAttribute('x', x - 5)
        marker.setAttribute('y', y - 5)
        marker.textContent = '(' + I + ' ' + J + ' ' + K + ')'
        }
    }

// Beim "`Loslassen"' des Atoms wird dessen Position auf das nächstgelegene Kohlenstoffatom gerundet. Das ist natürlich nur dann möglich, wenn sich dort nicht schon ein anderes Fremdatom befindet. Wenn von der Ausgangsposition aus schon ein "`zukünftiger"' Bewegungsverlauf gespeichert war, wird dieser überschrieben.

function drop() {
    var id = [I, J, K].join('-')

    if (! document.getElementById(id) || id == atom.id) {
        document.onmouseup = null

        if (id != atom.id)
            M.splice(++m, Number.MAX_VALUE, [atom.id, id])

        atom.setAttribute('cx', getX(I, J, K))
        atom.setAttribute('cy', getY(I, J, K))
        atom.id = id
        atom = null

        marker.setAttribute('x', -100)
        }
    }


// \codesection{Animierte Fremdatombewegung}

// Manuelle Bewegungen der Fremdatome oder solche, die vom implementierten Simulated-Annealing-Algorithmus durchgeführt wurden, können, wenn der Verlauf nicht anderweitig unterbrochen wurde, stets rückgängig gemacht oder wiederholt werden. Die einzelnen Transpositionen werden dabei von \verb|move()| animiert.

function move(id1, id2, n) {
    var object = document.getElementById(id1)
    object.id = id2

    clearInterval(step[id1])

    var x1 = parseFloat(object.getAttribute('cx'))
    var y1 = parseFloat(object.getAttribute('cy'))

    var ijk = id2.split('-').map(parseFloat)

    var x2 = getX.apply(this, ijk)
    var y2 = getY.apply(this, ijk)

    var i = 1
    var n = n || 40

    var dx = (x2 - x1) / n
    var dy = (y2 - y1) / n

    step[id2] = setInterval(function() {
        object.setAttribute('cx', x1 + i * dx)
        object.setAttribute('cy', y1 + i * dy)
        if (++i > n) clearInterval(step[id2])
        }, 10)
    }


// \codesection{Steuerung}

// Mit der Tastatur kann über \verb|type()| sowohl das Formular abgesendet (Enter), ein Schritt rückgängig gemacht ($\leftarrow$) sowie wiederholt werden ($\rightarrow$).

function type(event) {
    event = event || window.event

    switch (event.which || event.keyCode) {
        case 13: send(); break
        case 37: undo(); break
        case 39: redo(); break
        }
    }

// Letzgenannte Vorgänge werden durch \verb|undo()| und \verb|redo()| realisiert, die eine Animation starten und den Zeigerindex \verb|m| de- oder inkrementieren.

function undo(n) {
    if (m >= 0) {
        move(M[m][1], M[m][0], n)
        m--
        }
    }

function redo(n) {
    if (m < M.length - 1) {
        m++
        move(M[m][0], M[m][1], n)
        }
    }

// Eine automatische Wiederholung dieser Befehle, also das Abspielen eines Films, wird durch \verb|uundo()| und \verb|rredo()| veranlasst.

function uundo() {
    undo()
    if (m >= 0) play = setTimeout(uundo, 400)
    }

function rredo() {
    redo()
    if (m < M.length - 1) play = setTimeout(rredo, 400)
    }

// \verb|first()| und \verb|last()| bringen einen sofort zum einen oder anderen Ende der Bewegung.

function first() { while (m >= 0) undo(1) }

function last() { while (m < M.length - 1) redo(1) }

// Durch \verb|pause()| werden schließlich alle internen Prozesse gestoppt.

function pause() {
    window.stop()
    clearTimeout(play)
    clearTimeout(scan)
    form.anneal.value = ''
    form.scan.value = ''
    }
