;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Networks and Universal Mixing: Comparative compartmental disease models.
; (C) Christopher J Watts, 2020.
; Compares the universal mixing in a compartmental model of disease transmission
; with a social-network-based model.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

extensions [rnd csv]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

globals [
  Current-View ; Show compartments or network?

  last-seed-setup ; Random-number generator seeds
  last-seed-go

  sorted-people ; List, sorted by who. Avoids use of "ask people []", and hence generation of random permutations of people.

  susceptible ; Every model should have this disease state compartment.
  exposed ; First compartment after susceptible.
  compartments-sorted ; Used for updating in order, and one charts
  new-case-transition ; Define this to fit with your disease model.

  seed-infection-times ; List of time steps when seed infections occur.

  ; Disease Statistics
  net-cases-new
  net-cases-total
  net-cases-peak
  net-cases-peak-day
  net-infections-new
  net-infections-total
  net-infections-peak
  net-infections-peak-day

  um-cases-new
  um-cases-total
  um-cases-peak
  um-cases-peak-day
  um-infections-new
  um-infections-total
  um-infections-peak
  um-infections-peak-day

  ; Network statistics
  mean-degree
  median-degree
  min-degree
  max-degree
  var-degree
  assortativity
  network-density
  clustering-coefficient
  mean-cliquishness
  median-cliquishness
  min-cliquishness
  max-cliquishness
  mean-dos
  min-dos
  max-dos
  median-dos
  mean-closeness
  min-closeness
  max-closeness
  median-closeness
  net-diameter
  degree-centralization
  closeness-centralization

  num-components
  max-component-size


;  sim-time
  start-date-as-num
  end-date-as-num
]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

breed [people person]
undirected-link-breed [slinks slink]
breed [compartments compartment]
directed-link-breed [transitions transition]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

people-own [
  p-state ; Person's current disease state (compartment)
  p-maturity ; Time leave state, and which state person enters next.
  p-infection-time ; False if never infected. Otherwise, time they became infectious.
  p-circle ; Social circle
  p-cliquishness ; Local clustering coefficient.
  p-component ; Network component
  p-dos ; Average degree of separation from other nodes = 1 / Closeness
  p-tmp ; Used for DOS calc.
  p-reach
  p-closeness
]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

slinks-own [
  sl-contacts
]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

compartments-own [
  c-name
  c-color
  c-contents
  c-rel-inf ; relative infectiousness
]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

transitions-own [
  ts-weight
  ts-maturity-task
]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Setup

to setup
  clear-all
  ask patches [set pcolor white]

  set current-view view-of-compartments
  setup-disease-compartments
  setup-initial-plots
  change-view

  setup-rng "Seed-Setup"
  setup-population
  reset-sim
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-rng [given-variable-name]
  ifelse 0 = runresult given-variable-name [
    run (word "set last-" given-variable-name " " new-seed)
  ]
  [
    run (word "set last-" given-variable-name " " given-variable-name)
  ]
  random-seed runresult (word "last-" given-variable-name)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to switch-mechanism
  ifelse sim-mechanism = sim-mechanism-network [
    set sim-mechanism sim-mechanism-uni-mix
    reset-sim
  ]
  [
    set sim-mechanism sim-mechanism-network
    reset-sim
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report sim-mechanism-network
  report "Network"
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report sim-mechanism-uni-mix
  report "Universal Mixing"
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Disease transmission model

to setup-disease-compartments
  run (word "setup-disease-compartments-" compartments-design)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-disease-compartments-Susceptible-Infectious
  ; The simplest model: Susceptible-Infectious
  let mid-x mean (list max-pxcor min-pxcor)
  let x-step world-width / 8
  let y-step world-height / 4

  ; new-compartment name relative-infectiousness xcor ycor color
  let ds-S new-compartment "Susceptible" 0.0 mid-x (max-pycor - 1 * y-step) green
  let ds-I new-compartment "Infectious" 1.0 mid-x (max-pycor - 2 * y-step) red

  set susceptible ds-S
  set exposed ds-I
  set compartments-sorted (list ds-S ds-I)

  ; make-transition-to destination-compartment weight maturity-time-reporter
  ask ds-S [make-transition-to ds-I 1.0 [-> ]]

  ask ds-S [ask out-transition-to ds-I [set new-case-transition self]]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-disease-compartments-SIR
  ; The classic model introduced by Kermack & McKendrick (1927): Susceptible-Infectious-Removed
  let mid-x mean (list max-pxcor min-pxcor)
  let x-step world-width / 8
  let y-step world-height / 4

  ; new-compartment name relative-infectiousness xcor ycor color
  let ds-S new-compartment "Susceptible" 0.0 mid-x (max-pycor - 1 * y-step) green
  let ds-I new-compartment "Infectious" 1.0 mid-x (max-pycor - 2 * y-step) red
  let ds-R new-compartment "Removed" 0.0 (mid-x - 0 * x-step) (max-pycor - 3 * y-step) brown

  set susceptible ds-S
  set exposed ds-I
  set compartments-sorted (list ds-S ds-I ds-R)

  ; make-transition-to destination-compartment weight maturity-time-reporter
  ask ds-S [make-transition-to ds-I 1.0 [-> ]]
  ask ds-I [make-transition-to ds-R 1.0 [-> random-event-time-gamma 8.0 8.0]]

  ask ds-S [ask out-transition-to ds-I [set new-case-transition self]]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-disease-compartments-SEI3HRD
  ; Based on the Covid-19 model from the LSHTM (Davies et al. 2020)
  let mid-x mean (list max-pxcor min-pxcor)
  let x-step world-width / 8
  let y-step world-height / 8

  ; new-compartment name relative-infectiousness xcor ycor color
  ; new-compartment name relative-infectiousness xcor ycor
  let ds-S new-compartment "Susceptible" 0.0 mid-x (max-pycor - 1 * y-step) green
  let ds-E new-compartment "Exposed" 0.0 mid-x (max-pycor - 2 * y-step) yellow
  let ds-Ip new-compartment "I-Preclinical" 1.00 (mid-x + 1 * x-step) (max-pycor - 3 * y-step) orange
  let ds-Ic new-compartment "I-Clinical" 1.00 (mid-x + 1 * x-step) (max-pycor - 4 * y-step) red
  let ds-Is new-compartment "I-Subclinical" 0.50 (mid-x - 1 * x-step) (max-pycor - 4 * y-step) pink
  let ds-H new-compartment "Hospitalize?" 0.0 (mid-x + 1 * x-step) (max-pycor - 5 * y-step) sky
  let ds-ICU new-compartment "ICU" 0.0 (mid-x + 2 * x-step) (max-pycor - 6 * y-step) (violet - 1)
  let ds-Non-ICU new-compartment "Non-ICU" 0.0 (mid-x + 1 * x-step) (max-pycor - 6 * y-step) (violet + 2)
  let ds-R new-compartment "Recovered" 0.0 (mid-x - 1 * x-step) (max-pycor - 7 * y-step) brown
  let ds-D new-compartment "Dead" 0.0 (mid-x + 2 * x-step) (max-pycor - 7 * y-step) (gray - 4)

  set susceptible ds-S
  set exposed ds-E
  set compartments-sorted (list ds-S ds-E ds-Ip ds-Ic ds-Is ds-H ds-ICU ds-Non-ICU ds-R ds-D)

  ; make-transition-to destination-compartment weight maturity-time-reporter
  ask ds-S [make-transition-to ds-E 1.0 [-> ]]
  ask ds-E [make-transition-to ds-Ip (perc-symptomatic / 100) [-> random-event-time-gamma 4.0 4.0]]
  ask ds-E [make-transition-to ds-Is ((100 - perc-symptomatic) / 100) [-> random-event-time-gamma 4.0 4.0]]
  ask ds-Ip [make-transition-to ds-Ic 1.0 [-> random-event-time-gamma 1.5 4.0]]
  ask ds-Ic [make-transition-to ds-H 1.0 [-> random-event-time-gamma 3.5 3.5]]
  ask ds-Is [make-transition-to ds-R 1.0 [-> random-event-time-gamma 5.0 4.0]]

  ask ds-H [make-transition-to ds-R ((100 - perc-hospitalized) / 100) [-> 0]]
  ask ds-H [make-transition-to ds-ICU (perc-hospitalized * Perc-Need-ICU / 10000) [-> random-event-time-gamma 3.5 3.5]]
  ask ds-H [make-transition-to ds-Non-ICU (perc-hospitalized * (100 - Perc-Need-ICU) / 10000) [-> random-event-time-gamma 3.5 3.5]]
  ask ds-ICU [make-transition-to ds-R (Perc-Hospitalized-Recover / 100) [-> random-event-time-gamma 10.0 10.0]]
  ask ds-ICU [make-transition-to ds-D ((100 - Perc-Hospitalized-Recover) / 100) [-> random-event-time-gamma 11.5 11.5]]
  ask ds-Non-ICU [make-transition-to ds-R (Perc-Hospitalized-Recover / 100) [-> random-event-time-gamma 8.0 8.0]]
  ask ds-Non-ICU [make-transition-to ds-D ((100 - Perc-Hospitalized-Recover) / 100) [-> random-event-time-gamma 11.5 11.5]]

  ask ds-Ip [ask out-transition-to ds-Ic [set new-case-transition self]]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report new-compartment [given-name given-rel-inf given-x given-y given-color]
  let return-object nobody
  create-compartments 1 [
    set hidden? (Current-View != view-of-compartments)
    set return-object self
    set shape "square"
    set size 2
    set c-color given-color
    set color c-color
    setxy given-x given-y
    set label-color black
    set c-name given-name
    set c-rel-inf given-rel-inf
    set c-contents (turtle-set )
    relabel-compartment
  ]
  report return-object
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to relabel-compartment
  set label (word c-name ": " (count c-contents) "        ")
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report compartment-named [given-name]
  report min-one-of compartments with [c-name = given-name] [who]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to make-transition-to [given-state given-weight given-maturity-task]
  ; Run as a compartment
  create-transition-to given-state [
    set hidden? (Current-View != view-of-compartments)
    set color grey
    set thickness 0.2
    set ts-weight given-weight
    set ts-maturity-task given-maturity-task
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report random-event-time-gamma [mu k]
  ; Where mean=mu, shape=k, in days.
  ; Returns in days.
  report (random-gamma k (k / mu))
  report (int ((random-gamma k (k / mu)) / time-step)) * time-step
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-population
  setup-social-network
  foreach sorted-people [cur-person ->
    ask cur-person [
      set hidden? (current-view != view-of-network)
      set shape "person"
      set color red
      set size 1
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Social Networks

to setup-social-network
  run (word "setup-social-network-" network-architecture)
  calc-network-stats
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-slink [given-color]
  set hidden? (current-view != view-of-network)
  set color given-color
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-sorted-population
  set sorted-people []
  repeat population [
    create-people 1 [set sorted-people fput self sorted-people]
  ]
  set sorted-people reverse sorted-people
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-social-network-social-circles
  ; From Hamill & Gilbert (2009). http://jasss.soc.surrey.ac.uk/12/2/3.html
  setup-sorted-population
  foreach sorted-people [cur-person ->
    ask cur-person [
      setxy random-xcor random-ycor
      set p-circle 0
      set color red - 2
    ]
  ]
  foreach sorted-people [cur-person ->
    ask cur-person [
      create-slinks-with other (people in-radius soc-net-radius0) [setup-slink grey + 2]
    ]
  ]
  let group1 (turtle-set )
  ask n-of (perc-soc-net-group1 * (count people) / 100) people [
    set p-circle 1
    set group1 (turtle-set group1 self)
    set color red + 2
  ]
  foreach filter [p -> [1 = p-circle] of p] sorted-people [cur-person ->
    ask cur-person [
      create-slinks-with other (group1 in-radius soc-net-radius1) [setup-slink grey - 2]
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-social-network-erdos-renyi-random
  setup-sorted-population
  foreach sorted-people [ego ->
    ask ego [
      setxy random-xcor random-ycor
      foreach filter [p -> who > [who] of p] sorted-people [alter ->
        if er-density > random-float 1 [
          ask alter [
            create-slink-with ego [
              setup-slink grey + 2
            ]
          ]
        ]
      ]
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-social-network-ring
  setup-sorted-population
  layout-circle (sorted-people) (3 * (min (list world-width world-height)) / 8)
  if min-links-per-node = 0 [stop]
  let prev-nodes []
  let first-nodes []
  if min-links-per-node >= count people [user-message (word "Warning: \n\nmin-links-per-node (" min-links-per-node ") is >= Count People")]
  foreach sorted-people [ego ->
    if min-links-per-node > length first-nodes [
      set first-nodes fput ego first-nodes
    ]
    foreach prev-nodes [alter ->
      ask ego [
        create-slink-with alter [setup-slink grey + 2]
      ]
    ]
    set prev-nodes fput ego sublist prev-nodes 0 min (list (min-links-per-node - 1) length prev-nodes)
  ]
  foreach reverse first-nodes [ego ->
    foreach prev-nodes [alter ->
      ask ego [
        create-slink-with alter [setup-slink grey + 2]
      ]
    ]
    set prev-nodes fput ego sublist prev-nodes 0 min (list (min-links-per-node - 1) length prev-nodes)
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-social-network-strogatz-watts-small-world
  setup-social-network-ring
  let ego nobody
  let alter nobody
  foreach sort slinks [x ->
    if sw-rewire-chance > random-float 1.0 [
      set ego [end1] of x
      ask x [die]
      ask ego [
        create-slink-with one-of other people with [not slink-neighbor? myself] [setup-slink grey + 2]
      ]
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-social-network-barabasi-albert
;  nw:generate-preferential-attachment people slinks population min-links-per-node
;  layout-circle sort people (3 * (min (list world-width world-height)) / 8)
;  stop

  set sorted-people []

  repeat min-links-per-node [
    create-people 1 [
      foreach sorted-people [alter ->
        create-slink-with alter [setup-slink grey - 2]
      ]
      set sorted-people fput self sorted-people
    ]
  ]

;  let sorted-nodes reverse sort people
  repeat (population - min-links-per-node) [
    create-people 1 [
      foreach n-of min-links-per-node sorted-people [alter ->
        create-slink-with alter [
          setup-slink grey - 2
        ]
      ]
      set sorted-people fput self sorted-people
    ]
  ]
  layout-circle sorted-people (3 * (min (list world-width world-height)) / 8)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-grid-positions
  let cur-id 0
  let row-length sqrt population
  if row-length ^ 2 < population [set row-length row-length + 1]
  let x-step world-width / row-length
  let y-step world-height / row-length
  foreach sort people [p ->
    ask p [
      setxy (x-step * (cur-id mod row-length)) (y-step * (int (cur-id / row-length)))
      set cur-id cur-id + 1
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-social-network-2d-grid-neighbors [given-neighbors]
;  nw:generate-lattice-2d people slinks row-length row-length false
;  setup-grid-positions
;  stop

  set sorted-people []
  let row-length int sqrt population
  if row-length ^ 2 < population [set row-length row-length + 1]
  let cur-id 0
  let x-step world-width / (row-length + 1)
  let y-step world-height / (row-length + 1)
  foreach n-values row-length [? -> ?] [cur-row ->
    foreach n-values row-length [? -> ?] [cur-col ->
      create-people 1 [
        set sorted-people fput self sorted-people
;        show (list (x-step * (cur-col)) (y-step * (cur-row)))
        setxy (x-step * (cur-col)) (y-step * (cur-row))
      ]
    ]
  ]

  let neighbors-radius-pairs (list [4 1.01] [8 1.5] [12 2.01] [20 2.3] [24 2.85] [28 3.1] [36 3.2])

  let chosen-item position given-neighbors map [x -> first x] neighbors-radius-pairs
  if chosen-item = false [user-message "Fatal: Requested number of neighbors could not be found in list in setup-social-network-2d-grid-neighbors." stop]

  let chosen-radius x-step * last item chosen-item neighbors-radius-pairs
  foreach sorted-people [ego ->
    ask ego [
      create-slinks-with other people in-radius chosen-radius [
        setup-slink grey + 2
      ]
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-social-network-kissler-data
  if not file-exists? "Kissler_DataS1.csv" [user-message "FATAL: Cannot find file Kissler_DataS1.csv!" stop]
  if not file-exists? "Kissler_DataS2.csv" [user-message "FATAL: Cannot find file Kissler_DataS2.csv!" stop]
  let kd1 filter [r -> kd-max-distance >= item 3 r] csv:from-file "Kissler_DataS1.csv" ; Filter out all contacts at large distances.
  let node-ids sort remove-duplicates sentence (map [r -> item 1 r] kd1) (map [r -> item 2 r] kd1)
  set population length node-ids
  setup-sorted-population
  layout-circle (sorted-people) (7 * (min (list world-width world-height)) / 16)
  let who-adj min [who] of people
  foreach kd1 [r ->
    ask person (who-adj + position (item 1 r) node-ids) [
      if nobody = person (who-adj + position (item 2 r) node-ids) [print (word "Target node does not exist: " r)]
      ifelse slink-neighbor? person (who-adj + position (item 2 r) node-ids) [
        ask slink-with person (who-adj + position (item 2 r) node-ids) [set sl-contacts sl-contacts + 1]
      ]
      [
        create-slink-with person (who-adj + position (item 2 r) node-ids) [setup-slink lime set sl-contacts 1]
      ]
    ]
  ]

  let kd2 csv:from-file "Kissler_DataS2.csv"

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to change-view
  ; NB: Use of lists and sort. We don't want to use any random numbers here.
  ifelse current-view = view-of-network [
    set current-view view-of-compartments
  ]
  [
    if current-view = view-of-compartments [
      set current-view view-of-network
    ]
  ]
  if any? people [
    foreach sorted-people [p ->
      ask p [
        set hidden? (current-view != view-of-network)
        foreach sort (my-slinks with [end1 = myself]) [s ->
          ask s [
            set hidden? (current-view != view-of-network)
          ]
        ]
      ]
    ]
  ]
  foreach compartments-sorted [c ->
    ask c [
      set hidden? (Current-View != view-of-compartments)
      foreach sort (my-out-transitions) [t ->
        ask t [
          set hidden? (Current-View != view-of-compartments)
        ]
      ]
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report view-of-network
  report "Network"
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report view-of-compartments
  report "Compartments"
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Disease state transitions

to reset-sim
  reset-ticks
  set start-date-as-num parsed-date Start-Date
  set end-date-as-num parsed-date End-Date
  clear-relevant-plot

  setup-rng "Seed-Go"
  ; NB: If given the same RNG seed, this will initialize with the same seed infections after each reset.
  ; If you want different seed infections, but in the same network, then you need to generate and store them during Setup.
  initialize-population

  if sim-mechanism = sim-mechanism-network [
    set net-cases-new 0
    set net-cases-total 0
    set net-cases-peak 0
    set net-cases-peak-day 0
    set net-infections-new 0
    set net-infections-total 0
    set net-infections-peak 0
    set net-infections-peak-day 0
  ]

  if sim-mechanism = sim-mechanism-uni-mix [
    set um-cases-new 0
    set um-cases-total 0
    set um-cases-peak 0
    set um-cases-peak-day 0
    set um-infections-new 0
    set um-infections-total 0
    set um-infections-peak 0
    set um-infections-peak-day 0
  ]

  do-plots
  update-output
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to initialize-population
  ask people [
    set p-state susceptible
    set p-infection-time 0
    set p-maturity (list infinity)
    recolor-by-state
  ]
  foreach compartments-sorted [compa ->
    ask compa [set c-contents (turtle-set )]
  ]
  ask susceptible [set c-contents people]
  run schedule-seeds-n-d
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to schedule-seeds [num-seeds num-days]
  set seed-infection-times []
  foreach n-values num-days [d -> d] [cur-day ->
    foreach n-values num-seeds [n -> n] [cur-seed ->
      set seed-infection-times fput (cur-seed + round (cur-day / time-step)) seed-infection-times
    ]
  ]
  set seed-infection-times reverse seed-infection-times
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to become-infected
  let infected-person self
  ask susceptible [ask out-transition-to exposed [do-transition-by infected-person]]
  set p-infection-time sim-time
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to do-transition-by [given-person]
  ask end1 [set c-contents c-contents with [self != given-person] relabel-compartment]
  ask end2 [set c-contents (turtle-set c-contents given-person) relabel-compartment]
  if end1 = susceptible [if end2 = exposed [record-new-infection]]
  if self = new-case-transition [record-new-case]

  ask given-person [
    set p-state [end2] of myself
    recolor-by-state
    set p-maturity next-state-transition
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to record-new-infection
  if sim-mechanism = sim-mechanism-network [
    set net-infections-new net-infections-new + 1
    if net-infections-new > net-infections-peak [
      set net-infections-peak net-infections-new
      set net-infections-peak-day sim-time
    ]
    set net-infections-total net-infections-total + 1
    stop
  ]
  if sim-mechanism = sim-mechanism-uni-mix [
    set um-infections-new um-infections-new + 1
    if um-infections-new > um-infections-peak [
      set um-infections-peak um-infections-new
      set um-infections-peak-day sim-time
    ]
    set um-infections-total um-infections-total + 1
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to record-new-case
  if sim-mechanism = sim-mechanism-network [
    set net-cases-new net-cases-new + 1
    if net-cases-new > net-cases-peak [
      set net-cases-peak net-cases-new
      set net-cases-peak-day sim-time
    ]
    set net-cases-total net-cases-total + 1
    stop
  ]
  if sim-mechanism = sim-mechanism-uni-mix [
    set um-cases-new um-cases-new + 1
    if um-cases-new > um-cases-peak [
      set um-cases-peak um-cases-new
      set um-cases-peak-day sim-time
    ]
    set um-cases-total um-cases-total + 1
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report next-state-transition
  ; run by person on entering p-state
  if ([not any? my-out-transitions] of p-state) [report (list infinity)]
  let next-transition rnd:weighted-one-of ([my-out-transitions] of p-state) [ts-weight]
  report [(list (sim-time + runresult ts-maturity-task) self)] of next-transition
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to recolor-by-state
  set color [c-color] of p-state
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report contacts-per-time-step
  report time-step * contact-rate / 10
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to do-infections-in-network
  ; Susceptible ego contacts potentially infectious alter
  let inter-contact-time 1.0 / contacts-per-time-step
  ask [c-contents] of susceptible [
;    let num-contacts time-step * contact-rate / 10
;    set num-contacts (int num-contacts) + ifelse-value ((num-contacts mod 1) > random-float 1) [1] [0]
    ask slink-neighbors with [
      ifelse-value ([c-rel-inf > 0] of p-state) [time-step > random-exponential (inter-contact-time * (count slink-neighbors))] [false]
    ] [
;    ask n-of (random-poisson contacts-per-time-step) slink-neighbors [
      if [susceptible != p-state] of myself [stop] ; Ego may have already become infectious, in which case don't waste computer time.
      if [c-rel-inf > 0] of p-state [ ; Alter is infectious.
        if [c-rel-inf * susceptibility > random-float 100] of p-state [
          ask myself [become-infected]
        ]
      ]
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to do-infections-by-universal-mixing
  let inf-chance 1.0 - exp (0 - force-of-infection * time-step)
  ask [c-contents with [inf-chance > random-float 1]] of susceptible [ ; No random-binomial!
    become-infected
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report force-of-infection
  report susceptibility * contacts-per-time-step * ((sum [c-rel-inf * count c-contents] of compartments) / count people) / 100
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to do-maturation
  foreach but-first compartments-sorted [compa ->
    ask [c-contents with [sim-time >= first p-maturity]] of compa [
      ask last p-maturity [do-transition-by myself]
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report infinity
  ; Arbitrarily large number.
  report 2147483648 ; = 2 ^ 31
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to go
  if sim-time = end-date-as-num - start-date-as-num [stop]

  if 0 = sim-time mod 1 [ ; New day.
    if sim-mechanism = sim-mechanism-network [set net-cases-new 0 set net-infections-new 0]
    if sim-mechanism = sim-mechanism-uni-mix [set um-cases-new 0 set um-infections-new 0]
  ]

  ; Do seed infections (if any scheduled)
  while [ifelse-value (not empty? seed-infection-times) [sim-time = first seed-infection-times] [false]] [
    ask susceptible [
      if any? c-contents [
        ask one-of c-contents [become-infected]
      ]
    ]
    set seed-infection-times but-first seed-infection-times
  ]

  ; Do infections.
  if sim-mechanism = sim-mechanism-uni-mix [do-infections-by-universal-mixing]
  if sim-mechanism = sim-mechanism-network [do-infections-in-network]

  ; Do other state transitions in order.
  do-maturation

  tick
  do-plots
  update-output
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report sim-time
  report ticks * time-step
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Network statistics

to calc-network-stats
  calc-components
  if not any? slinks [stop]
  set mean-degree mean [count my-slinks] of people
  set median-degree median [count my-slinks] of people
  set min-degree min [count my-slinks] of people
  set max-degree max [count my-slinks] of people
  set var-degree variance [count my-slinks] of people
  set degree-centralization (count people) * (max-degree - mean-degree) / (((count people) - 1) * ((count people) - 2))

  set network-density 200 * (count slinks) / ((count people) * (-1 + count people))

  if Calculate-Slow-Metrics? [
    calc-cliquishness
    ;if num-components = 1 [calc-dos]
    calc-dos
  ]

  let assort-list [(list ([count my-slinks] of end1 ) ([count my-slinks] of end2 ))] of slinks
;  let assort-list [(list (count my-slinks) (mean [count my-slinks] of slink-neighbors))] of people
  let x-bar mean map [p -> first p] assort-list
  let y-bar mean map [p -> last p] assort-list
  let s_x sqrt mean map [p -> ((first p) - x-bar) ^ 2] assort-list
  let s_y sqrt mean map [p -> ((last p) - y-bar) ^ 2] assort-list
  ifelse (s_x * s_y) > 0 [
    set assortativity (mean map [p -> ((first p) - x-bar) * ((last p) - y-bar) ] assort-list) / (s_x * s_y)
  ]
  [
    set assortativity 1.0
  ]

  do-network-plots
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to calc-cliquishness
  let num-triplets 0
  let num-triangles 0

  foreach sorted-people [ego ->
    ask ego [
      ifelse 2 < count my-slinks [
        set p-cliquishness 200 * (sum [count (my-slinks with [end1 = myself]) with [[slink-neighbor? ego] of end2] ] of slink-neighbors) /
        ((count my-slinks) * (-1 + count my-slinks))
      ]
      [
        set p-cliquishness 0
      ]

      ; Count triplets and triangles
      foreach sort (slink-neighbors with [who > [who] of myself]) [alter ->
        foreach sort (slink-neighbors with [who > [who] of alter]) [tertius ->
          set num-triplets num-triplets + 1
          if [slink-neighbor? tertius] of alter [
            set num-triangles num-triangles + 1
          ]
        ]
      ]
    ]
  ]
  ifelse num-triplets = 0 [
    set clustering-coefficient 0
  ]
  [
    set clustering-coefficient 100 * num-triangles / num-triplets
  ]

  set mean-cliquishness mean [p-cliquishness] of people
  set median-cliquishness median [p-cliquishness] of people
  set min-cliquishness min [p-cliquishness] of people
  set max-cliquishness max [p-cliquishness] of people
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to calc-components
  let node-stack []
  let ego nobody
  set num-components 0
  set max-component-size 0
  let cur-size 0

  foreach sorted-people [ ?1 ->
    ask ?1 [set p-component 0]
  ]

  foreach sorted-people [ ?1 ->
    if [p-component = 0] of ?1 [
      set cur-size 0
      set num-components 1 + num-components
      set node-stack fput ?1 node-stack
    ]

    while [not empty? node-stack] [
      set ego first node-stack
      set node-stack but-first node-stack
      if [p-component = 0] of ego [
        ask ego [
          set cur-size cur-size + 1
          if cur-size > max-component-size [set max-component-size cur-size]
          set p-component num-components
          ask slink-neighbors [
            if p-component = 0 [
              set node-stack fput self node-stack
            ]
          ]
        ]
      ]
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to calc-dos
  ; From Ulrik Brandes's (2001) betweenness algorithm
  ; NB: Meaningless if network not fully connected. (i.e. some distances between nodes would be infinite!)
  let Q [] ; Queue (FIFO)
  let v nobody
  let maxdos 0
  set net-diameter -1

  foreach sorted-people [ego ->
    ask ego [
      foreach sorted-people [alter ->
        ask alter [set p-tmp -1]
      ]
      set p-tmp 0
      set Q []
      set Q lput self Q
      while [length Q > 0] [
        set v first Q
        set Q but-first Q
        ask [slink-neighbors] of v [
          if p-tmp < 0 [
            set Q lput self Q
            set p-tmp (1 + [p-tmp] of v)
          ]
        ]
      ]

      set p-reach count people with [p-tmp >= 0]
      set p-dos (sum [p-tmp] of people) / p-reach
      set p-closeness ifelse-value (p-dos <= 0) [-1] [(p-reach - 1) / (p-dos * p-reach)]
      set maxdos max [p-tmp] of people
      if (maxdos > net-diameter) [ set net-diameter maxdos ]
    ]
  ]

  set mean-dos mean [p-dos] of people
  set median-dos median [p-dos] of people
  set min-dos min [p-dos] of people
  set max-dos max [p-dos] of people

  set mean-closeness mean [p-closeness] of people
  set median-closeness median [p-closeness] of people
  set min-closeness min [p-closeness] of people
  set max-closeness max [p-closeness] of people

  set closeness-centralization (count people) * (max-closeness - mean-closeness) * (2 * (count people) - 3) / (((count people) - 1) * ((count people) - 2))
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plots and outputs

to do-network-plots
  set-current-plot "Degree of Connectivity"
  set-plot-x-range 0 (1 + max [count my-slinks] of people)
  set-plot-pen-interval 1
  histogram [count my-slinks] of people
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-initial-plots
  set-current-plot "Network"
  setup-extra-plot-pens

  set-current-plot "Universal Mixing"
  setup-extra-plot-pens
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-extra-plot-pens
  clear-plot
  let pen-name "Susceptible"
;  set-current-plot-pen pen-name
;  set-plot-pen-color [c-color] of susceptible
  foreach compartments-sorted [cur-comp ->
    set pen-name [c-name] of cur-comp
    create-temporary-plot-pen pen-name
    set-current-plot-pen pen-name
    set-plot-pen-color [c-color] of cur-comp
    set-plot-pen-mode 0
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to clear-relevant-plot
  if sim-mechanism = sim-mechanism-network [
    set-current-plot "Network"
    setup-extra-plot-pens
    set-current-plot "New Infections"
    set-current-plot-pen "Network"
    plot-pen-reset
    set-current-plot "New Cases"
    set-current-plot-pen "Network"
    plot-pen-reset
    stop
  ]
  if sim-mechanism = sim-mechanism-uni-mix [
    set-current-plot "Universal Mixing"
    setup-extra-plot-pens
    set-current-plot "New Infections"
    set-current-plot-pen "Universal Mix."
    plot-pen-reset
    set-current-plot "New Cases"
    set-current-plot-pen "Universal Mix."
    plot-pen-reset
    stop
  ]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to do-plots
  if sim-mechanism = sim-mechanism-network [
    set-current-plot "Network"
    do-compartments-plot

    set-current-plot "New Infections"
    set-current-plot-pen "Network"
    plotxy sim-time (as-perc net-infections-new)

    set-current-plot "New Cases"
    set-current-plot-pen "Network"
    plotxy sim-time (as-perc net-cases-new)

  ]
  if sim-mechanism = sim-mechanism-uni-mix [
    set-current-plot "Universal Mixing"
    do-compartments-plot

    set-current-plot "New Infections"
    set-current-plot-pen "Universal Mix."
    plotxy sim-time (as-perc um-infections-new)

    set-current-plot "New Cases"
    set-current-plot-pen "Universal Mix."
    plotxy sim-time (as-perc um-cases-new)
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to do-compartments-plot
  let pen-name ""
  foreach compartments-sorted [cur-comp ->
    set pen-name [c-name] of cur-comp
    set-current-plot-pen pen-name
    plotxy sim-time [as-perc count c-contents] of cur-comp
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to update-output
  clear-output
  if Update-Output? [
    output-print (word "Day = " sim-time ":")
    foreach compartments-sorted [compa ->
      output-print [(word c-name " : " (count c-contents) " (" (precision (as-perc count c-contents) 1) "%)")] of compa
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report as-perc [given-count]
  report 100 * given-count / count people
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report parsed-date [given-string]
  ; Convert date string in format yyyy-m-d
  ; to a number. Used with start-date for parsing school holiday dates.
  let sep1 position "-" given-string
  let sep2 sep1 + 1 + position "-" substring given-string (sep1 + 1) (length given-string)
  let y cnum substring given-string 0 sep1
  let m cnum substring given-string (sep1 + 1) sep2
  let d cnum substring given-string (sep2 + 1) (length given-string)
  ; Treat 2020 as year 0. 2020 was a leapyear.
  let ret sum n-values (y - 2020) [x -> ifelse-value (0 = x mod 4) [366] [365]]
  set ret ret + ifelse-value (0 = y mod 4) [sum sublist [31 29 31 30 31 30 31 31 30 31 30 31] 0 (m - 1)] [sum sublist [31 28 31 30 31 30 31 31 30 31 30 31] 0 (m - 1)]
  ; 2020-1-1 is day 0
  set ret ret + d - 1
  report ret ; Return in days.
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report date-string [given-sim-time]
  ; Given-sim-time is in days.
  let tmp-sum start-date-as-num + int (given-sim-time )
  ; Day0 is 2020-1-1
  let cur-year 2020
  set tmp-sum tmp-sum - ifelse-value (0 = cur-year mod 4) [366] [365]
  while [tmp-sum > 0] [
    set cur-year cur-year + 1
    set tmp-sum tmp-sum - ifelse-value (0 = cur-year mod 4) [366] [365]
  ]
  set tmp-sum tmp-sum + ifelse-value (0 = cur-year mod 4) [366] [365]
  let cur-month 0
  set tmp-sum tmp-sum - item cur-month [31 29 31 30 31 30 31 31 30 31 30 31]
  while [tmp-sum > 0] [
    set cur-month cur-month + 1
    set tmp-sum tmp-sum - item cur-month [31 29 31 30 31 30 31 31 30 31 30 31]
  ]
  set tmp-sum tmp-sum + item cur-month [31 29 31 30 31 30 31 31 30 31 30 31]
  report (word cur-year "-" (cur-month + 1) "-" (1 + int tmp-sum)) ; yyyy-m-d
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report cnum [given-string]
  let ret 0
  let cur-item 0
  repeat length given-string [
    set ret 10 * ret + position (item cur-item given-string) ["0" "1" "2" "3" "4" "5" "6" "7" "8" "9"]
    set cur-item cur-item + 1
  ]
  report ret
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@#$#@#$#@
GRAPHICS-WINDOW
465
10
857
403
-1
-1
12.0
1
10
1
1
1
0
0
0
1
0
31
0
31
1
1
1
ticks
30.0

INPUTBOX
10
170
110
230
Population
392.0
1
0
Number

INPUTBOX
10
660
125
720
Soc-Net-Radius0
2.1
1
0
Number

INPUTBOX
130
660
245
720
Soc-Net-Radius1
3.5
1
0
Number

SLIDER
245
660
427
693
Perc-Soc-Net-Group1
Perc-Soc-Net-Group1
0
100
25.0
1
1
%
HORIZONTAL

CHOOSER
280
420
440
465
Schedule-Seeds-n-d
Schedule-Seeds-n-d
"schedule-seeds 1 1" "schedule-seeds 1 7"
1

SLIDER
10
420
260
453
Contact-Rate
Contact-Rate
0
500
250.0
1
1
Contacts Per 10 Days
HORIZONTAL

MONITOR
10
455
117
500
Contacts Per Day
contact-rate / 10
17
1
11

INPUTBOX
115
170
210
230
Time-Step
0.25
1
0
Number

MONITOR
120
455
262
500
Contacts Per Time-Step
contacts-per-time-step
17
1
11

SLIDER
205
505
405
538
Perc-Symptomatic
Perc-Symptomatic
0
100
33.0
1
1
%
HORIZONTAL

SLIDER
10
540
200
573
Perc-Hospitalized
Perc-Hospitalized
0
100
10.0
1
1
%
HORIZONTAL

BUTTON
220
185
322
218
Change View
change-view
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
330
180
460
225
NIL
Current-View
17
1
11

PLOT
870
10
1170
210
Network
Day
% of Pop
0.0
1.0
0.0
100.0
true
true
"" ""
PENS
"Susceptible" 1.0 0 -10899396 true "" ""

PLOT
870
215
1170
415
Universal Mixing
Day
% of Pop
0.0
1.0
0.0
100.0
true
true
"" ""
PENS
"Susceptible" 1.0 0 -10899396 true "" ""

BUTTON
220
70
290
103
NIL
Setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
315
285
385
318
NIL
Go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
390
285
460
318
Go x 1
Go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
10
10
430
60
Networks Versus Universal Mixing: Comparative Models of Disease Transmission
20
0.0
1

MONITOR
245
330
350
375
Day
Sim-Time
17
1
11

MONITOR
870
420
1170
465
NIL
Seed-Infection-Times
17
1
11

MONITOR
870
465
987
510
NIL
Force-Of-Infection
17
1
11

OUTPUT
1175
350
1415
595
11

CHOOSER
345
230
460
275
Sim-Mechanism
Sim-Mechanism
"Network" "Universal Mixing"
0

SLIDER
205
575
405
608
Perc-Hospitalized-Recover
Perc-Hospitalized-Recover
0
100
50.0
1
1
%
HORIZONTAL

SLIDER
10
505
200
538
Susceptibility
Susceptibility
0
100
20.0
1
1
%
HORIZONTAL

BUTTON
220
285
290
318
Reset
reset-sim
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
220
240
340
273
Switch Mechanism
switch-mechanism
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
1175
600
1317
633
Update-Output?
Update-Output?
1
1
-1000

MONITOR
300
70
460
115
Population Density (%)
100 * population / count patches
2
1
11

MONITOR
300
120
460
165
Network Density (%)
Network-Density
1
1
11

MONITOR
465
525
522
570
Mean
mean-degree
1
1
11

MONITOR
585
525
642
570
Median
median-degree
17
1
11

MONITOR
525
525
582
570
Min
min-degree
17
1
11

MONITOR
645
525
702
570
Max
Max-degree
17
1
11

TEXTBOX
465
505
615
523
Degree of connectivity:
13
0.0
1

MONITOR
705
525
762
570
Variance
var-degree
1
1
11

MONITOR
1375
50
1440
95
Total
Net-Cases-Total
17
1
11

TEXTBOX
1175
30
1250
48
Infections:
13
0.0
1

TEXTBOX
1315
30
1385
48
Cases:
13
0.0
1

MONITOR
1315
50
1375
95
New
Net-Cases-New
17
1
11

MONITOR
1315
95
1375
140
Peak
net-cases-peak
17
1
11

MONITOR
1375
95
1440
140
Peak Day
net-cases-peak-day
17
1
11

MONITOR
1175
50
1235
95
New
net-infections-new
17
1
11

MONITOR
1235
50
1300
95
Total
net-infections-total
17
1
11

MONITOR
1175
95
1235
140
Peak
net-infections-peak
17
1
11

MONITOR
1235
95
1300
140
Peak Day
net-infections-peak-day
17
1
11

MONITOR
1175
235
1235
280
New
um-infections-new
17
1
11

MONITOR
1235
235
1300
280
Total
um-infections-total
17
1
11

MONITOR
1175
280
1235
325
Peak
um-infections-peak
17
1
11

MONITOR
1235
280
1300
325
Peak Day
um-infections-peak-day
17
1
11

MONITOR
1315
235
1375
280
New
um-cases-new
17
1
11

MONITOR
1375
235
1440
280
Total
um-cases-total
17
1
11

MONITOR
1315
280
1375
325
Peak
um-cases-peak
17
1
11

MONITOR
1375
280
1440
325
Peak Day
um-cases-peak-day
17
1
11

TEXTBOX
1265
10
1330
28
Network:
13
0.0
1

TEXTBOX
1255
195
1405
213
Universal Mixing:
13
0.0
1

TEXTBOX
1175
215
1255
233
Infections:
13
0.0
1

TEXTBOX
1315
215
1465
233
Cases:
13
0.0
1

MONITOR
465
575
550
620
NIL
Assortativity
3
1
11

CHOOSER
10
120
212
165
Network-Architecture
Network-Architecture
"Social-Circles" "Erdos-Renyi-Random" "Barabasi-Albert" "Ring" "Strogatz-Watts-Small-World" "2D-Grid-Neighbors 4" "2D-Grid-Neighbors 8" "2D-Grid-Neighbors 12" "2D-Grid-Neighbors 20" "2D-Grid-Neighbors 28" "2D-Grid-Neighbors 36" "Kissler-Data"
0

INPUTBOX
145
775
260
835
ER-Density
0.014
1
0
Number

INPUTBOX
280
775
395
835
Min-Links-Per-Node
3.0
1
0
Number

MONITOR
555
430
637
475
NIL
Count SLinks
17
1
11

MONITOR
640
430
737
475
SLinks per Node
(count slinks) / (count people)
3
1
11

MONITOR
465
430
552
475
NIL
Count People
17
1
11

INPUTBOX
10
775
125
835
SW-Rewire-Chance
0.01
1
0
Number

TEXTBOX
15
390
350
415
Disease Transmission Parameters:
16
0.0
1

TEXTBOX
10
615
160
633
Network Parameters:
16
0.0
1

TEXTBOX
10
640
160
658
Social Circles:
13
0.0
1

TEXTBOX
145
738
265
773
Erdos-Renyi Random:
13
0.0
1

TEXTBOX
10
735
120
775
Strogatz-Watts Small World:
13
0.0
1

TEXTBOX
280
735
400
765
Various (Ring, SW, BA):
13
0.0
1

PLOT
555
575
760
725
Degree of Connectivity
Neighbors per Node
Count Nodes
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" ""

INPUTBOX
10
235
110
295
Start-Date
2020-1-29
1
0
String

INPUTBOX
115
235
210
295
End-Date
2020-05-29
1
0
String

MONITOR
355
330
460
375
Date (Y-M-D)
date-string sim-time
17
1
11

CHOOSER
10
70
210
115
Compartments-Design
Compartments-Design
"Susceptible-Infectious" "SIR" "SEI3HRD"
2

TEXTBOX
1180
145
1440
186
NB: (New) infections and (new) cases are evemts, changes in stock levels.
11
0.0
1

INPUTBOX
10
970
162
1030
Seed-Setup
0.0
1
0
Number

INPUTBOX
10
1035
162
1095
Seed-Go
0.0
1
0
Number

BUTTON
165
995
227
1028
Clear
set seed-setup 0
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
230
995
337
1028
Use Last Seed
set seed-setup last-seed-setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
165
1060
227
1093
Clear
set seed-go 0
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
230
1060
337
1093
Use Last Seed
set seed-go last-seed-go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
10
940
240
976
Random-Number Generation:
16
0.0
1

TEXTBOX
10
1100
345
1205
Remember: \nDifferent network structures may make very different use of random numbers.\nUniversal mixing will make different use of random numbers to transmission in a network.\nSo using the common seed numbers may not be meaningful.
11
0.0
1

MONITOR
90
330
207
375
Run-Length (Days)
end-date-as-num - start-date-as-num
1
1
11

PLOT
870
525
1170
675
New Infections
Day
% of Pop
0.0
1.0
0.0
1.0
true
true
"" ""
PENS
"Network" 1.0 0 -11221820 true "" ""
"Universal Mix." 1.0 0 -5825686 true "" ""

PLOT
870
675
1170
825
New Cases
Day
% of Pop
0.0
1.0
0.0
1.0
true
true
"" ""
PENS
"Network" 1.0 0 -11221820 true "" ""
"Universal Mix." 1.0 0 -5825686 true "" ""

MONITOR
465
620
550
665
Density (%)
Network-Density
2
1
11

BUTTON
625
740
687
773
NIL
Setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
460
810
517
855
Mean
Mean-cliquishness
1
1
11

MONITOR
520
810
577
855
Min
Min-cliquishness
1
1
11

MONITOR
640
810
697
855
Max
max-cliquishness
1
1
11

MONITOR
580
810
637
855
Median
Median-cliquishness
1
1
11

MONITOR
460
730
617
775
Clustering-Coefficient (%)
Clustering-Coefficient
1
1
11

TEXTBOX
460
785
735
803
Cliquishness (Local Clustering) (% of Triplets):
13
0.0
1

TEXTBOX
465
485
615
503
Network Statistics:
16
0.0
1

SLIDER
205
540
405
573
Perc-Need-ICU
Perc-Need-ICU
0
100
33.0
1
1
%
HORIZONTAL

TEXTBOX
460
860
675
878
Network Components:
13
0.0
1

MONITOR
460
880
567
925
NIL
Num-Components
17
1
11

MONITOR
570
880
697
925
NIL
Max-Component-Size
17
1
11

MONITOR
460
980
517
1025
Mean
mean-dos
1
1
11

TEXTBOX
460
940
610
958
Degrees of Separation:
13
0.0
1

MONITOR
520
980
577
1025
Min
min-dos
1
1
11

MONITOR
580
980
637
1025
Median
median-dos
1
1
11

MONITOR
640
980
697
1025
Max
max-dos
1
1
11

MONITOR
460
1030
547
1075
NIL
Net-Diameter
17
1
11

MONITOR
550
1030
697
1075
NIL
Closeness-Centralization
3
1
11

MONITOR
555
1080
687
1125
NIL
Degree-Centralization
3
1
11

SWITCH
625
490
807
523
Calculate-Slow-Metrics?
Calculate-Slow-Metrics?
0
1
-1000

BUTTON
690
740
827
773
NIL
Calc-Network-Stats
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
10
865
125
925
KD-Max-Distance
2.0
1
0
Number

TEXTBOX
10
845
115
863
Kissler-Data:
13
0.0
1

TEXTBOX
460
960
650
980
(Meaningless if Num-Components > 1)
11
0.0
1

BUTTON
315
380
377
413
Spring
repeat 10 [\nlayout-spring people slinks 1.0 1.0 1.0\n]
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

@#$#@#$#@
# NETWORKS AND UNIVERSAL MIXING

## WHAT IS IT?

This program simulates the transmission of a disease using one of two mechanisms. In both mechanisms, people susceptible to the disease become infected via social contact with those who are infectious. One mechanism assumes those susceptible to the disease tend to mix uniformly with those infectious (the "universal mixing" assumption). The second mechanism assumes that social contacts are constrained by social or spatial networks. Spatially and socially remote people cannot directly infect each other. The program is designed to allow comparisons between universal mixing and various social network structures.

This program (C) Christopher J. Watts, 2020.

### Compartmental models in Epidemiology

Compartmental models of disease transmission are a common tool in epidemiology. The members of a population are distributed among a number of compartments, representing the states of a disease. Infections occur when members of one compartment ("Susceptible") make social contacts with members of another ("Infectious"). As a result of becoming infected, people leave the "Susceptible" compartment and flow into another, eventually arriving in the "Infectious" compartment. Most models, including the best known one, the S-I-R model of Kermack & McKendrick (1927), also include flow out of the "Infectious" compartment at some fixed rate over time. It is also common to place compartments between Susceptible and Infectious, e.g. representing an incubation period, when a person has been exposed to the disease, but is not yet infectious. E.g. in the SEIR model, compartments represent being "Susceptible" to the disease, being "Exposed", or infected but not yet infectious, then being "Infectious", and finally "Removed" from being able to contact and infect other people. This last state may include those who have self-isolated (quarantined), those who have been isolated in hospital, those who have recovered and are now immune, and those who have died. More complex models may wish to distinguish between these states of removal, e.g. for estimating how many people will need hospital beds, or how many people will die. Other variations include multiple compartments for Infectious states, e.g. for Covid-19 many infectious people display relatively few or no symptoms. We may wish to distinguish their degree of infectiousness and also their behaviour from those who do show symptoms.

### Social networks and the assumption of universal mixing

As a simplification, most compartmental models assume universal mixing between susceptible and infectious people, i.e. anyone's contacts can be with anyone else. If your model simulates an entire country (e.g. the UK), susceptibles in the far north of the country can be infected by people in the far south. This can seem implausible.

In reality, people live in a geographical space, and epidemics reflect this. As well as following spatial structures, people's contacts reflect their social networks. Clearly our contacts are constrained to the spatial and social networks we live in. Exactly what structure our social networks have is an open question. In addition, it may be unclear which types of social relations are most relevant to a particular disease. E.g. different social relations are involved in a sexually transmitted disease to a pneumonic disease that spreads through coughing or sneezing. Hence, we may need to try multiple types of network structure to make robust forecasts of the likely spread and impact of a disease.

This program allows comparisons between transmission via universal mixing and transmission via a social network. A variety of social network structures ("architectures") are provided, though others could be added (e.g. from the nw extension). In addition, different disease models (SI, SIR, SEIR etc.), with various combinations of compartments, can be tested.

## HOW IT WORKS

A number of compartments are created, representing disease states, and linked by directed links representing state transitions. A social network of people is created. Each person is assigned a disease state. A schedule is created by which, at specified time points, a number of randomly chosen people will become infected (from sources outside the model). Each time step, social contacts between any susceptible and infectious people are simulated. For the "Network" mechanism, the choice of recipients of these contacts is restricted to the susceptible person's social network neighbors. For the "Universal Mixing" mechanism, a function is used to calculate the "force of infection" on a susceptible individual. The force of infection is the rate at which infection events occur for that susceptible person as a result of social contacts with any infectious people in the population. 

### What is common to "Networks" and "Universal Mixing"

The program is designed so that the two mechanisms operate on the same size population, with the same seed infections, the same attributes of the disease (susceptibility, symptomatic rate, relative infectiousness etc.), and the same contact rate. Both simulation mechanisms include stochastic processes.

### Where the mechanisms differ

While the two transmission mechanisms share the same contact rate, they differ in the question of whom those contacts can be with. 

In "Universal Mixing", any infectious person can contribute to the chance of a particular susceptible person becoming infected in the current time step. Because any of them can contribute, there is no need to identify them or simulate individual interactions. The "force of infection" procedure just needs to count how many people there are in each infectious state. This saves on the compuutational work required for the simulation, and this mechanism should run faster than with "Network". In "Networks", more stochastic processes are explicitly simulated, requiring more random number generation, and hence more computer time. 

This has the consequence that, if that person becomes infected we can say who exactly infected them, and study the attributes of those who infect more people than most ("superspreaders").

The user may feel that the results under "Universal Mixing" are less meaningful. Under universal mixing, an infectious person contributes to a susceptible person's force of infection for every time step the two remain in their respective states. It is as if the susceptible could meet anyone in the population at any time. In reality, our interaction opportunities are more constrained than this. The day you met someone in the far north cannot be the same day you met someone in the south.

In "Network", we sample contacts from the susceptible person's social link neighbours. If no one among their neighbours is currently infectious, then the chance of being infected in the current time step is zero, no matter how prevalent the disease is elsewhere in the network. In addition, the outcomes of interactions in adjacent time steps will tend to be highly correlated, especially if the social network is highly clustered. If none of your neighbours are infectious, then probably most of your neighbours' neighbours are uninfectious, and so your neighbours are likely to remain uninfectious for some more time steps. On the other hand, once the disease reaches your cluster, you are likely to join them in an infected state soon after. 



## HOW TO USE IT

Choose a compartmental model design (e.g. "SIR"), and a network architecture (e.g. "Social Circles"). Click "Setup". A population will be created, given social network links, and distributed across the NetLogo World (e.g. given random coordinates, or organised into a ring shape). At the same time, a disease compartmental model is created with disease states represented by a breed of turtles ("compartments"), and state transitions represented by directed links. Initially these are hidden. You can switch to viewing them (and making the people in their network hidden), by clicking on the button "Change View". Choose which Sim-Mechanism you want to run first ("Network" or "Universal Mixing"). Then click "Go". When you have finished, you can either "Reset" the population, and repeat the simulation using the same mechanism, or click the "Switch Mechanism" button. This preserves the results from the mechanism run with so far, but resets the population. Now, clicking "Go" will run using the alternative mechanism, and update its corresponding plot and monitors. You can then compare the plots and statistics for the two mechanisms.

## THINGS TO NOTICE

Compared to a largely spatial network (e.g. "Social Circles"), universal mixing tends to produce sharper, taller peaks in infections and cases. However, the Erdos-Renyi random network is quite similar to universal mixing, despite it having much less than 100% network density. (A complete network, with 100% of the possible pairs of nodes having a link present, is the equivalent of universal mixing.)

## THINGS TO TRY

How do epidemic outcomes (total cases, peak cases, time to peak cases, proportion of the population still susceptible at the end) vary with network attributes, such as density, degree distribution, cluster coefficient, and assortativity? How do outcomes for individual nodes (infected / not infected, time of infection) vary with the nodes' attributes (degree of connectivity, cliquishness)?

Where network architectures vary in their epidemic outcomes, how sensitive is this behaviour to the structure of the disease model? Would a choice of different compartments affect this, or different tranisition probabilities (weights) and times? 

## EXTENDING THE MODEL

### More network architectures!

The nw extension should allow you to quickly add more architectures.  (Though see the note below about performance for larger networks.)

### More compartmental model designs!

The code for the compartmental models is intended to be easy to understand and copy. To make your own compartmental model, copy one of the existing "setup-disease-compartments-" procedures. Add lines for a new-compartment and make-new-transition. Add your new compartments to the "compartments-sorted" list. Make there the global variable "exposed" is always set to a compartment. Infection events are defined as the transition from Susceptible to whatever state Exposed is defined as. (In the S-I model, the Exposed variable is set to the "Infectious" compartment.) You may need to alter the definition of Cases. E.g. in an SEIR model, to the transition from Exposed to Infectious. In the SEI3R model, this is the transition from Infectious-Preclinical to Infectious-Clinical, and thus excludes those who are asymptomatic (Infectious-Subclinical). 

### More agent attributes

Age, sex, health conditions, occupation, hobbies. E.g. what if older people are more likely to show symptoms, require hospitalization, or die?

### Death and other behaviour-changing events

Notice how, in the force of infection function, the chance of a contact being with an infectious person is assumed to be a multiple of the infectious proportion of people __in the entire population__. Is it really plausible that we make social contact with the dead? Or make contacts at the same rate as before with the symptomatic, the hospitalized, the newly recovered (and probably not yet fit enough to go to work / school)? We spend time on social contacts, and if there are now fewer people to interact with (because some are isolating, dead, etc.), that leaves more time to interact with the remaining people, some of whom are infectious. How could we reflect this potential for reallocation of time and intensification of contacts on fewer others?

## NETLOGO FEATURES

For some network architectures, I have used my own code, and not the nw extension. This was partly due to the nw versions being too slow for larger networks (e.g. >10000 nodes).

Use is made of the rnd extension for weighted sampling. (Thanks again to Nicolas Payette!)


## RELATED MODELS

Since the Coivd-19 pandemic forced modellers to stay home, there have probably emerged more mathematical and computer simulation models than will ever be countable. 

This program owes much to an attempt to replicate in a NetLogo agent-based model the compartmental model of Covid-19 written by Adam Kucharski and colleagues at the London School of Hygiene and Tropical Medicines (Davies et al. 2020a, 2020b). Kucharski's program is a stochastic compartmental model, written in R and C++, with age-based symptomatic rates.

All compartmental models are descended from the SIR model of Kermack & McKendrick (1927).

Disease epidemic models are only one type of model of a spreading phenomenon. The diffusion of innovations is often modelled using some of the same techniques (rightly or wrongly). See Watts & Gilbert (2014, especially Ch.s 2, 3).

## CREDITS AND REFERENCES

### Epidemic models

Davies, N. G., Kucharski, A. J., Eggo, R. M., Gimma, A., Edmunds, W. J., Jombart, T., . . . Liu, Y. (2020a). Effects of non-pharmaceutical interventions on COVID-19 cases, deaths, and demand  for  hospital  services  in  the  UK:  a  modelling  study.  The  Lancet  Public  Health,  5(7),  e375-e385. doi:10.1016/S2468-2667(20)30133-X 

Davies, N. G., Kucharski, A. J., Eggo, R. M., Gimma, A., Edmunds, W. J., Jombart, T., . . . Liu, Y. (2020b). Supplement to Davies et al.(2020a). https://www.thelancet.com/cms/10.1016/S2468-2667(20)30133-X/attachment/cee85e76-cffb-42e5-97b6-06a7e1e2379a/mmc1.pdf

Kermack, William Ogilvy, & A. G. McKendrick (1927) "A contribution to the mathematical theory of epidemics". Procedings of the Royal Society A. Volume 115, Issue 772. https://doi.org/10.1098/rspa.1927.0118

Watts, Christopher and Gilbert, Nigel (2014) “Simulating Innovation: Computer- based Tools for Rethinking Innovation”. Cheltenham, UK and Northampton, MA, USA: Edward Elgar.

### Networks

Barabási, A.-L., & Albert, R. (1999). “Emergence of scaling in random networks.” Science, 286(5439), 509-512. doi: 10.1126/science.286.5439.509

Erdős, P., & Rényi, A. (1959). “On random graphs I.” Publicationes Mathematicae, 6, 290-297.

Hamill, L., & Gilbert, N. (2009). “Social Circles: A Simple Structure for Agent-Based Social Network Models.” Jasss-the Journal of Artificial Societies and Social Simulation, 12(2). doi: 3

Klepac, Petra, Stephen Kissler, & Julia Gog (2018) "Contagion! The BBC Four Pandemic – The model behind the documentary". Epidemics Volume 24, September 2018, Pages 49-59.

https://github.com/skissler/haslemere

Watts, D. J., & Strogatz, S. H. (1998). “Collective dynamics of ‘small-world’ networks.” Nature, 393(6684), 440-442. doi: 10.1038/30918

## OUR TERMS AND CONDITIONS OF USE

If you use this program in your work, please cite the author, Christopher J. Watts, and the URL for the website from which you downloaded the program and the date on which you downloaded it. 

This program is free software: you can redistribute it and/or modify it under the terms of version 3 of the GNU General Public License as published by the Free Software Foundation. See http://www.gnu.org/licenses/ for more details. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

Those wishing to make commercial use of the program should contact the program’s author, Christopher J. Watts, to discuss terms.

This program has been designed to be run using NetLogo version 6.1.1, which can be obtained from http://ccl.northwestern.edu/netlogo/
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment-UM" repetitions="50" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>timer</metric>
    <metric>count people</metric>
    <metric>um-infections-new</metric>
    <metric>um-infections-total</metric>
    <metric>um-infections-peak</metric>
    <metric>um-infections-peak-day</metric>
    <metric>um-cases-new</metric>
    <metric>um-cases-total</metric>
    <metric>um-cases-peak</metric>
    <metric>um-cases-peak-day</metric>
    <metric>[count c-contents] of susceptible</metric>
    <metric>[as-perc count c-contents] of susceptible</metric>
    <enumeratedValueSet variable="Compartments-Design">
      <value value="&quot;Susceptible-Infectious&quot;"/>
      <value value="&quot;SIR&quot;"/>
      <value value="&quot;SEI3HRD&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Architecture">
      <value value="&quot;Social-Circles&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Population">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Sim-Mechanism">
      <value value="&quot;Universal Mixing&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Calculate-Slow-Metrics?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Start-Date">
      <value value="&quot;2020-1-29&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="End-Date">
      <value value="&quot;2021-1-29&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Schedule-Seeds-n-d">
      <value value="&quot;schedule-seeds 1 7&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Contact-Rate">
      <value value="50"/>
      <value value="75"/>
      <value value="100"/>
      <value value="125"/>
      <value value="150"/>
      <value value="175"/>
      <value value="200"/>
      <value value="225"/>
      <value value="250"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Susceptibility">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Symptomatic">
      <value value="33"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Hospitalized">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Hospitalized-Recover">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Soc-Net-Radius0">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Soc-Net-Radius1">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Soc-Net-Group1">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Min-Links-Per-Node">
      <value value="14"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SW-Rewire-Chance">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ER-Density">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Time-Step">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Update-Output?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Seed-Setup">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Seed-Go">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-Net" repetitions="50" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>timer</metric>
    <metric>count people</metric>
    <metric>net-infections-new</metric>
    <metric>net-infections-total</metric>
    <metric>net-infections-peak</metric>
    <metric>net-infections-peak-day</metric>
    <metric>net-cases-new</metric>
    <metric>net-cases-total</metric>
    <metric>net-cases-peak</metric>
    <metric>net-cases-peak-day</metric>
    <metric>[count c-contents] of susceptible</metric>
    <metric>[as-perc count c-contents] of susceptible</metric>
    <metric>count slinks</metric>
    <metric>(count slinks) / count people</metric>
    <metric>mean-degree</metric>
    <metric>median-degree</metric>
    <metric>min-degree</metric>
    <metric>max-degree</metric>
    <metric>var-degree</metric>
    <metric>assortativity</metric>
    <metric>network-density</metric>
    <metric>clustering-coefficient</metric>
    <metric>mean-cliquishness</metric>
    <metric>median-cliquishness</metric>
    <metric>min-cliquishness</metric>
    <metric>max-cliquishness</metric>
    <metric>mean-dos</metric>
    <metric>min-dos</metric>
    <metric>max-dos</metric>
    <metric>median-dos</metric>
    <metric>mean-closeness</metric>
    <metric>min-closeness</metric>
    <metric>max-closeness</metric>
    <metric>median-closeness</metric>
    <metric>net-diameter</metric>
    <metric>degree-centralization</metric>
    <metric>closeness-centralization</metric>
    <metric>num-components</metric>
    <metric>max-component-size</metric>
    <enumeratedValueSet variable="Compartments-Design">
      <value value="&quot;Susceptible-Infectious&quot;"/>
      <value value="&quot;SIR&quot;"/>
      <value value="&quot;SEI3HRD&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Architecture">
      <value value="&quot;Social-Circles&quot;"/>
      <value value="&quot;Erdos-Renyi-Random&quot;"/>
      <value value="&quot;Barabasi-Albert&quot;"/>
      <value value="&quot;Ring&quot;"/>
      <value value="&quot;Strogatz-Watts-Small-World&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Population">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Sim-Mechanism">
      <value value="&quot;Network&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Calculate-Slow-Metrics?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Start-Date">
      <value value="&quot;2020-1-29&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="End-Date">
      <value value="&quot;2020-10-31&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Schedule-Seeds-n-d">
      <value value="&quot;schedule-seeds 1 7&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Contact-Rate">
      <value value="50"/>
      <value value="75"/>
      <value value="100"/>
      <value value="125"/>
      <value value="150"/>
      <value value="175"/>
      <value value="200"/>
      <value value="225"/>
      <value value="250"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Susceptibility">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Symptomatic">
      <value value="33"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Hospitalized">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Hospitalized-Recover">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Soc-Net-Radius0">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Soc-Net-Radius1">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Soc-Net-Group1">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Min-Links-Per-Node">
      <value value="14"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SW-Rewire-Chance">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ER-Density">
      <value value="0.028"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Time-Step">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Update-Output?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Seed-Setup">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Seed-Go">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-2D-Grid" repetitions="50" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>timer</metric>
    <metric>count people</metric>
    <metric>net-infections-new</metric>
    <metric>net-infections-total</metric>
    <metric>net-infections-peak</metric>
    <metric>net-infections-peak-day</metric>
    <metric>net-cases-new</metric>
    <metric>net-cases-total</metric>
    <metric>net-cases-peak</metric>
    <metric>net-cases-peak-day</metric>
    <metric>[count c-contents] of susceptible</metric>
    <metric>[as-perc count c-contents] of susceptible</metric>
    <metric>count slinks</metric>
    <metric>(count slinks) / count people</metric>
    <metric>mean-degree</metric>
    <metric>median-degree</metric>
    <metric>min-degree</metric>
    <metric>max-degree</metric>
    <metric>var-degree</metric>
    <metric>assortativity</metric>
    <metric>network-density</metric>
    <metric>clustering-coefficient</metric>
    <metric>mean-cliquishness</metric>
    <metric>median-cliquishness</metric>
    <metric>min-cliquishness</metric>
    <metric>max-cliquishness</metric>
    <metric>mean-dos</metric>
    <metric>min-dos</metric>
    <metric>max-dos</metric>
    <metric>median-dos</metric>
    <metric>mean-closeness</metric>
    <metric>min-closeness</metric>
    <metric>max-closeness</metric>
    <metric>median-closeness</metric>
    <metric>net-diameter</metric>
    <metric>degree-centralization</metric>
    <metric>closeness-centralization</metric>
    <metric>num-components</metric>
    <metric>max-component-size</metric>
    <enumeratedValueSet variable="Compartments-Design">
      <value value="&quot;Susceptible-Infectious&quot;"/>
      <value value="&quot;SIR&quot;"/>
      <value value="&quot;SEI3HRD&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Architecture">
      <value value="&quot;2D-Grid-Neighbors 4&quot;"/>
      <value value="&quot;2D-Grid-Neighbors 8&quot;"/>
      <value value="&quot;2D-Grid-Neighbors 12&quot;"/>
      <value value="&quot;2D-Grid-Neighbors 20&quot;"/>
      <value value="&quot;2D-Grid-Neighbors 28&quot;"/>
      <value value="&quot;2D-Grid-Neighbors 36&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Population">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Sim-Mechanism">
      <value value="&quot;Network&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Calculate-Slow-Metrics?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Start-Date">
      <value value="&quot;2020-1-29&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="End-Date">
      <value value="&quot;2021-1-29&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Schedule-Seeds-n-d">
      <value value="&quot;schedule-seeds 1 7&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Contact-Rate">
      <value value="50"/>
      <value value="75"/>
      <value value="100"/>
      <value value="125"/>
      <value value="150"/>
      <value value="175"/>
      <value value="200"/>
      <value value="225"/>
      <value value="250"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Susceptibility">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Symptomatic">
      <value value="33"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Hospitalized">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Hospitalized-Recover">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Soc-Net-Radius0">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Soc-Net-Radius1">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Soc-Net-Group1">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Min-Links-Per-Node">
      <value value="14"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SW-Rewire-Chance">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ER-Density">
      <value value="0.028"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Time-Step">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Update-Output?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Seed-Setup">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Seed-Go">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-UM-392" repetitions="50" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>timer</metric>
    <metric>count people</metric>
    <metric>um-infections-new</metric>
    <metric>um-infections-total</metric>
    <metric>um-infections-peak</metric>
    <metric>um-infections-peak-day</metric>
    <metric>um-cases-new</metric>
    <metric>um-cases-total</metric>
    <metric>um-cases-peak</metric>
    <metric>um-cases-peak-day</metric>
    <metric>[count c-contents] of susceptible</metric>
    <metric>[as-perc count c-contents] of susceptible</metric>
    <enumeratedValueSet variable="Compartments-Design">
      <value value="&quot;Susceptible-Infectious&quot;"/>
      <value value="&quot;SIR&quot;"/>
      <value value="&quot;SEI3HRD&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Architecture">
      <value value="&quot;Social-Circles&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Population">
      <value value="392"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Sim-Mechanism">
      <value value="&quot;Universal Mixing&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Calculate-Slow-Metrics?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Start-Date">
      <value value="&quot;2020-1-29&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="End-Date">
      <value value="&quot;2021-1-29&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Schedule-Seeds-n-d">
      <value value="&quot;schedule-seeds 1 7&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Contact-Rate">
      <value value="50"/>
      <value value="75"/>
      <value value="100"/>
      <value value="125"/>
      <value value="150"/>
      <value value="175"/>
      <value value="200"/>
      <value value="225"/>
      <value value="250"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Susceptibility">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Symptomatic">
      <value value="33"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Hospitalized">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Hospitalized-Recover">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Soc-Net-Radius0">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Soc-Net-Radius1">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Soc-Net-Group1">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Min-Links-Per-Node">
      <value value="14"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SW-Rewire-Chance">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ER-Density">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Time-Step">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Update-Output?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Seed-Setup">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Seed-Go">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-Net-392" repetitions="50" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>timer</metric>
    <metric>count people</metric>
    <metric>net-infections-new</metric>
    <metric>net-infections-total</metric>
    <metric>net-infections-peak</metric>
    <metric>net-infections-peak-day</metric>
    <metric>net-cases-new</metric>
    <metric>net-cases-total</metric>
    <metric>net-cases-peak</metric>
    <metric>net-cases-peak-day</metric>
    <metric>[count c-contents] of susceptible</metric>
    <metric>[as-perc count c-contents] of susceptible</metric>
    <metric>count slinks</metric>
    <metric>(count slinks) / count people</metric>
    <metric>mean-degree</metric>
    <metric>median-degree</metric>
    <metric>min-degree</metric>
    <metric>max-degree</metric>
    <metric>var-degree</metric>
    <metric>assortativity</metric>
    <metric>network-density</metric>
    <metric>clustering-coefficient</metric>
    <metric>mean-cliquishness</metric>
    <metric>median-cliquishness</metric>
    <metric>min-cliquishness</metric>
    <metric>max-cliquishness</metric>
    <metric>mean-dos</metric>
    <metric>min-dos</metric>
    <metric>max-dos</metric>
    <metric>median-dos</metric>
    <metric>mean-closeness</metric>
    <metric>min-closeness</metric>
    <metric>max-closeness</metric>
    <metric>median-closeness</metric>
    <metric>net-diameter</metric>
    <metric>degree-centralization</metric>
    <metric>closeness-centralization</metric>
    <metric>num-components</metric>
    <metric>max-component-size</metric>
    <enumeratedValueSet variable="Compartments-Design">
      <value value="&quot;Susceptible-Infectious&quot;"/>
      <value value="&quot;SIR&quot;"/>
      <value value="&quot;SEI3HRD&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Network-Architecture">
      <value value="&quot;2D-Grid-Neighbors 4&quot;"/>
      <value value="&quot;2D-Grid-Neighbors 8&quot;"/>
      <value value="&quot;Kissler-Data&quot;"/>
      <value value="&quot;Social-Circles&quot;"/>
      <value value="&quot;Erdos-Renyi-Random&quot;"/>
      <value value="&quot;Barabasi-Albert&quot;"/>
      <value value="&quot;Ring&quot;"/>
      <value value="&quot;Strogatz-Watts-Small-World&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Population">
      <value value="392"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Sim-Mechanism">
      <value value="&quot;Network&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Calculate-Slow-Metrics?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Start-Date">
      <value value="&quot;2020-1-29&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="End-Date">
      <value value="&quot;2021-1-29&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Schedule-Seeds-n-d">
      <value value="&quot;schedule-seeds 1 7&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Contact-Rate">
      <value value="50"/>
      <value value="75"/>
      <value value="100"/>
      <value value="125"/>
      <value value="150"/>
      <value value="175"/>
      <value value="200"/>
      <value value="225"/>
      <value value="250"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Susceptibility">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Symptomatic">
      <value value="33"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Hospitalized">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Hospitalized-Recover">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Soc-Net-Radius0">
      <value value="2.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Soc-Net-Radius1">
      <value value="3.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Perc-Soc-Net-Group1">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Min-Links-Per-Node">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SW-Rewire-Chance">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ER-Density">
      <value value="0.014"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Time-Step">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Update-Output?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Seed-Setup">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Seed-Go">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
1
@#$#@#$#@
