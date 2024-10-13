//
//  ContentView.swift
//  FluidSimulator
//
//  Created by Nick Venanzi on 9/29/24.
//

import SwiftUI

struct ContentView: View {
    
    var engine = PhysicsEngine()
    @State var count = 0
    
    var body: some View {
        Button {
            engine.iterate()
            count += 1
        } label: {
            Text("Iterate: \(count)")
        }

    }
}

#Preview {
    ContentView()
}
